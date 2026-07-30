// Parallel variant of uniid_fold: fold one shard and write its accumulator and
// compressed proof to disk. Run one instance per shard (e.g. a job array).
//
// Place in Nova-Scotia/examples/ and run:
//   cargo run --release --example uniid_fold_shard -- <r1cs> <wasm> <shard.json> <out_prefix>
//
// Writes <out_prefix>.acc, <out_prefix>.proof, <out_prefix>.meta.

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
    if args.len() < 5 {
        eprintln!("usage: uniid_fold_shard <r1cs> <wasm> <shard.json> <out_prefix>");
        std::process::exit(1);
    }
    let r1cs_path = PathBuf::from(&args[1]);
    let wasm_path = PathBuf::from(&args[2]);
    let out_prefix = &args[4];

    let v: Value = serde_json::from_str(&fs::read_to_string(&args[3]).unwrap()).unwrap();
    let root_s = v["root"].as_str().unwrap().to_string();
    let acc_iv_s = v["acc_iv"].as_str().unwrap().to_string();
    let expected_s = v["expected_acc"].as_str().unwrap().to_string();
    let shard = v["shard"].as_str().unwrap_or("shard").to_string();
    let steps = v["steps"].as_array().unwrap();
    let n = steps.len();
    println!("shard {}  steps {}", shard, n);

    let mut private_inputs = Vec::with_capacity(n);
    for s in steps {
        let mut m = HashMap::new();
        for key in ["chrom", "pos", "refCode", "altCode", "pathElements", "pathIndices"] {
            m.insert(key.to_string(), s[key].clone());
        }
        private_inputs.push(m);
    }

    let start_public_input = [fe(&acc_iv_s), fe(&root_s)];
    let z0_secondary = [F::<G2>::from(0)];
    let r1cs = load_r1cs::<G1, G2>(&FileLocation::PathBuf(r1cs_path));
    let pp: PublicParams<G1, G2, _, _> = create_public_params(r1cs.clone());

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

    let (z, _) = recursive_snark
        .verify(&pp, n, &start_public_input, &z0_secondary)
        .expect("recursive verify failed");
    let acc_final = z[0];
    assert!(acc_final == fe(&expected_s), "accumulator != expected_acc");
    assert!(z[1] == fe(&root_s), "root not preserved");

    let (pk, vk) = CompressedSNARK::<_, _, _, _, S<G1>, S<G2>>::setup(&pp).unwrap();
    let t2 = Instant::now();
    let compressed =
        CompressedSNARK::<_, _, _, _, S<G1>, S<G2>>::prove(&pp, &pk, &recursive_snark).unwrap();
    let cprove_s = t2.elapsed().as_secs_f64();
    assert!(compressed
        .verify(&vk, n, start_public_input.to_vec(), z0_secondary.to_vec())
        .is_ok());

    fs::write(format!("{}.acc", out_prefix), format!("{:?}\n", acc_final)).unwrap();
    fs::write(format!("{}.proof", out_prefix), bincode::serialize(&compressed).unwrap()).unwrap();
    fs::write(
        format!("{}.meta", out_prefix),
        format!("shard\t{}\nnSteps\t{}\nroot\t{}\nacc\t{:?}\nfold_s\t{:.1}\ncompress_prove_s\t{:.1}\n",
            shard, n, root_s, acc_final, fold_s, cprove_s),
    )
    .unwrap();
    println!("fold {:.1} s; wrote {}.acc/.proof/.meta", fold_s, out_prefix);
}
