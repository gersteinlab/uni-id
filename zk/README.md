# uni-id zero-knowledge proofs

Zero-knowledge proofs over the LP (common-variant) layer of the uni-id partition.

- **Single-variant membership**: prove one variant belongs to the allowed common-variant set, revealing nothing about which variant.
- **All-variants proof**: prove that every variant in an LP set belongs to the allowed set and that the set folds to a published commitment, in one succinct proof, revealing no variant.

Both prove membership against a Poseidon Merkle tree against a given allowed set of variants. Hashing in circuits is Poseidon; the single-variant proof uses Groth16 (circom + snarkjs), the all-variants proof uses Nova folding (via Nova-Scotia).

## Layout

```
circuits/   mtree.circom         shared Merkle templates
            membership.circom    single-variant proof (Groth16)
            lpfold.circom        fold step for the all-variants proof (Nova)
scripts/    build_membership_input.js   tree + one-variant witness
            build_fold_inputs.js        tree + all paths + running accumulator
            build_shard_inputs.js       parallel: per-shard inputs
            combine_clp.js              parallel: fold shard accumulators -> commitment
rust/       uniid_fold.rs               all-variants fold driver
            uniid_fold_shard.rs         parallel: fold one shard
            combine_clp.rs              parallel: verify all shard proofs
```

## Requirements

- circom 2.x, Node 18+, and the JS deps: `npm install` (installs circomlibjs, circomlib, snarkjs).
- For the all-variants proof: Rust and [Nova-Scotia](https://github.com/nalinbhardwaj/Nova-Scotia). Nova-Scotia pulls arkworks crates that need a recent Rust (1.90+). The `rust/*.rs` files are Nova-Scotia examples; see below.

## Input

`allowed.bed.gz`: bgzip BED of the allowed set, columns `chrom start end ref alt` (0-based start). Variants are referenced by key `chrom:pos:ref:alt` with 1-based `pos = start + 1` (chr prefix optional). Leaf encoding: `Poseidon(chromInt, pos, refCode, altCode)`, where alleles are base-6 encoded (A=1, C=2, G=3, T=4, N=5). The tree pads to a power of two with `Poseidon(0,0,0,0)`.

`depth` in the circuits must match the tree (23 for a ~6.3M-leaf set). Edit the `main` line in `membership.circom` / `lpfold.circom` for a different set size.

## Single-variant membership

```bash
npm install
cd circuits
circom membership.circom --r1cs --wasm --sym -l ../node_modules/circomlib/circuits
cd ..

# witness for one variant (omit the key to prove the middle leaf)
node scripts/build_membership_input.js allowed.bed.gz all . 22:10562724:T:C

# trusted setup (2^15 covers depth 23; reuse across proofs)
npx snarkjs powersoftau new bn128 15 pot0.ptau
npx snarkjs powersoftau contribute pot0.ptau pot1.ptau --name=t1 -e="random"
npx snarkjs powersoftau prepare phase2 pot1.ptau pot.ptau
npx snarkjs groth16 setup circuits/membership.r1cs pot.ptau mem0.zkey
npx snarkjs zkey contribute mem0.zkey mem.zkey --name=k -e="random"
npx snarkjs zkey export verificationkey mem.zkey vkey.json

# prove and verify
node circuits/membership_js/generate_witness.js circuits/membership_js/membership.wasm input.json w.wtns
npx snarkjs groth16 prove mem.zkey w.wtns proof.json public.json
npx snarkjs groth16 verify vkey.json public.json proof.json
```

`public.json` contains only the root. Verification succeeds only if the variant is in the tree.

## All-variants proof

Set up Nova-Scotia once:

```bash
git clone https://github.com/nalinbhardwaj/Nova-Scotia.git
cp rust/uniid_fold.rs Nova-Scotia/examples/
# add to Nova-Scotia/Cargo.toml [dependencies]: ff = "0.13", serde_json = "1", bincode = "1.3"
```

Compile the fold step and build inputs from a key list (one `chrom:pos:ref:alt` per line, in fold order):

```bash
cd circuits
circom lpfold.circom --r1cs --wasm --sym -l ../node_modules/circomlib/circuits
cd ..
node scripts/build_fold_inputs.js allowed.bed.gz all . lp_keys.txt   # -> fold_inputs.json
```

Fold, verify, compress:

```bash
cd Nova-Scotia
cargo run --release --example uniid_fold -- \
  ../circuits/lpfold.r1cs ../circuits/lpfold_js/lpfold.wasm ../fold_inputs.json
```

## Whole-genome (parallel)

A single fold chain over millions of variants is sequential and slow. Divide the ordered key list into fixed-size shards and fold each independently from the same seed; the commitment is a single Poseidon fold over the shard accumulators. Trust rests only on the final commitment, not on the intermediate accumulators.

```bash
# 1. split the ordered key list into shards
mkdir shards
split -d -a 4 -l 50000 --additional-suffix=.keys lp_keys.txt shards/shard_

# 2. build per-shard inputs (one tree build; writes shards/shard_*.json and shard_acc.tsv)
node scripts/build_shard_inputs.js allowed.bed.gz shards shards
#    also write a manifest: shard_id \t nSteps
for f in shards/shard_*.keys; do
  id=$(basename "$f" .keys | sed 's/shard_//'); printf "%s\t%d\n" "$id" "$(wc -l < "$f")"
done > manifest.tsv

# 3. fold each shard (run in parallel, e.g. a job array over the shard ids)
cp rust/uniid_fold_shard.rs rust/combine_clp.rs Nova-Scotia/examples/
mkdir proofs
cargo run --release --example uniid_fold_shard -- \
  circuits/lpfold.r1cs circuits/lpfold_js/lpfold.wasm shards/shard_0000.json proofs/shard_0000

# 4. verify all shard proofs and compute the commitment
cargo run --release --example combine_clp -- \
  circuits/lpfold.r1cs circuits/lpfold_js/lpfold.wasm shards proofs manifest.tsv
node scripts/combine_clp.js shard_acc.tsv    # prints C_LP
```
