// Parallel: verify every shard's compressed proof. The commitment C_LP is folded
// from the shard accumulators by combine_clp.js.
//
// Place in Nova-Scotia/examples/ and run:
//   cargo run --release --example combine_clp -- <r1cs> <wasm> <shards_dir> <proofs_dir> <manifest.tsv>
//
// manifest.tsv rows: shard_id \t nSteps

use std::{fs, path::PathBuf, time::Instant};

use ff::PrimeField;
use nova_scotia::{circom::reader::load_r1cs, create_public_params, FileLocation, F, S};
use nova_snark::{provider, CompressedSNARK, PublicParams};

type G1 = provider::bn256_grumpkin::bn256::Point;
type G2 = provider::bn256_grumpkin::grumpkin::Point;

fn fe(s: &str) -> F<G1> {
    F::<G1>::from_str_vartime(s).expect("bad field element")
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 6 {
        eprintln!("usage: combine_clp <r1cs> <wasm> <shards_dir> <proofs_dir> <manifest.tsv>");
        std::process::exit(1);
    }
    let r1cs_path = PathBuf::from(&a[1]);
    let shards_dir = PathBuf::from(&a[3]);
    let proofs_dir = PathBuf::from(&a[4]);

    let mut ids: Vec<(String, usize)> = fs::read_to_string(&a[5])
        .unwrap()
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let mut it = l.split('\t');
            (it.next().unwrap().to_string(), it.next().unwrap().trim().parse().unwrap())
        })
        .collect();
    ids.sort_by(|x, y| x.0.cmp(&y.0));

    let r1cs = load_r1cs::<G1, G2>(&FileLocation::PathBuf(r1cs_path));
    let pp: PublicParams<G1, G2, _, _> = create_public_params(r1cs);
    let (_pk, vk) = CompressedSNARK::<_, _, _, _, S<G1>, S<G2>>::setup(&pp).unwrap();
    let z0_secondary = vec![F::<G2>::from(0)];

    let mut all_ok = true;
    let t = Instant::now();
    for (id, n) in &ids {
        let v: serde_json::Value = serde_json::from_str(
            &fs::read_to_string(shards_dir.join(format!("shard_{}.json", id))).unwrap(),
        )
        .unwrap();
        let start = vec![fe(v["acc_iv"].as_str().unwrap()), fe(v["root"].as_str().unwrap())];
        let bytes = fs::read(proofs_dir.join(format!("shard_{}.proof", id))).unwrap();
        let compressed: CompressedSNARK<G1, G2, _, _, S<G1>, S<G2>> =
            bincode::deserialize(&bytes).unwrap();
        let ok = compressed.verify(&vk, *n, start, z0_secondary.clone()).is_ok();
        println!("shard {}  steps {}  {}", id, n, if ok { "OK" } else { "FAIL" });
        all_ok &= ok;
    }
    println!("\n{} shard proofs verified in {:?}: {}", ids.len(), t.elapsed(), all_ok);
    assert!(all_ok);
}
