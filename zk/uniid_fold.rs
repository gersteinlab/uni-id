// All-variants proof: fold the LP set into one recursive proof (Nova, BN254/Grumpkin),
// verify, check the accumulator equals the committed value, and compress.
//
// Place in Nova-Scotia/examples/ and run:
//   cargo run --release --example uniid_fold -- <r1cs> <wasm> <fold_inputs.json>
//
// z_0 = [acc_iv, root]. Each step checks membership against root and folds one
// variant. Proof exists only if every variant is in the tree and the final
// accumulator equals fold_inputs.expected_acc.

use std::{collections::HashMap, fs, path::PathBuf, time::Instant};

use ff::PrimeField;
use nova_scotia::{
    circom::reader::load_r1cs, create_public_params, create_recursive_circuit, FileLocation, F, S,
};
use nova_snark::{provider, CompressedSNARK, PublicParams};
use serde_json::Value;

type G1 = provider::bn256_grumpkin::bn256::Point;
type G2 = provider::bn256_grumpkin::grumpkin::Point;

fn fe(s: &str) -> F<G1> {
    F::<G1>::from_str_vartime(s).expect("bad field element")
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        eprintln!("usage: uniid_fold <r1cs> <wasm> <fold_inputs.json>");
        std::process::exit(1);
    }
    let r1cs_path = PathBuf::from(&args[1]);
    let wasm_path = PathBuf::from(&args[2]);

    let v: Value = serde_json::from_str(&fs::read_to_string(&args[3]).unwrap()).unwrap();
    let root_s = v["root"].as_str().unwrap();
    let acc_iv_s = v["acc_iv"].as_str().unwrap();
    let expected_s = v["expected_acc"].as_str().unwrap();
    let steps = v["steps"].as_array().unwrap();
    let n = steps.len();
    println!("variants {}  depth {}", n, v["depth"]);

    let mut private_inputs = Vec::with_capacity(n);
    for s in steps {
        let mut m = HashMap::new();
        for key in ["chrom", "pos", "refCode", "altCode", "pathElements", "pathIndices"] {
            m.insert(key.to_string(), s[key].clone());
        }
        private_inputs.push(m);
    }

    let start_public_input = [fe(acc_iv_s), fe(root_s)];
    let z0_secondary = [F::<G2>::from(0)];
    let r1cs = load_r1cs::<G1, G2>(&FileLocation::PathBuf(r1cs_path));
    let pp: PublicParams<G1, G2, _, _> = create_public_params(r1cs.clone());
    println!("constraints/step: primary {} secondary {}", pp.num_constraints().0, pp.num_constraints().1);

    let t = Instant::now();
    let recursive_snark = create_recursive_circuit(
        FileLocation::PathBuf(wasm_path),
        r1cs.clone(),
        private_inputs,
        start_public_input.to_vec(),
        &pp,
    )
    .unwrap();
    let fold_s = t.elapsed().as_secs_f64();
    println!("fold {:.1} s ({:.4} s/step)", fold_s, fold_s / n as f64);

    let t = Instant::now();
    let (z, _) = recursive_snark
        .verify(&pp, n, &start_public_input, &z0_secondary)
        .expect("recursive verify failed");
    println!("recursive verify {:.3} s", t.elapsed().as_secs_f64());
    assert!(z[0] == fe(expected_s), "accumulator != expected_acc");
    assert!(z[1] == fe(root_s), "root not preserved");
    println!("accumulator matches committed value; root preserved");

    let (pk, vk) = CompressedSNARK::<_, _, _, _, S<G1>, S<G2>>::setup(&pp).unwrap();
    let t = Instant::now();
    let compressed =
        CompressedSNARK::<_, _, _, _, S<G1>, S<G2>>::prove(&pp, &pk, &recursive_snark).unwrap();
    println!("compress prove {:.1} s", t.elapsed().as_secs_f64());

    let t = Instant::now();
    let ok = compressed
        .verify(&vk, n, start_public_input.to_vec(), z0_secondary.to_vec())
        .is_ok();
    println!("compressed verify {:.3} s: {}", t.elapsed().as_secs_f64(), ok);
    assert!(ok);
}
