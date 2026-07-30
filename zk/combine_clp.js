// Fold the per-shard accumulators into the single published commitment C_LP,
// using the same Poseidon as the circuits.
//
//   node combine_clp.js <shard_acc.tsv>
//
// shard_acc.tsv rows: shard_id \t nSteps \t acc   (from build_shard_inputs.js)
// C_LP = Poseidon fold over acc values in shard-id order, seeded at 0.

const fs = require("fs");
const { buildPoseidon } = require("circomlibjs");

const SEED = 0n;

async function main() {
  const tsv = process.argv[2];
  if (!tsv) {
    console.error("usage: node combine_clp.js <shard_acc.tsv>");
    process.exit(1);
  }
  const rows = fs.readFileSync(tsv, "utf8")
    .split("\n").map((l) => l.trim()).filter((l) => l.length)
    .map((l) => l.split("\t"))
    .sort((a, b) => a[0].localeCompare(b[0]));

  const poseidon = await buildPoseidon();
  const F = poseidon.F;
  const H = (arr) => F.toObject(poseidon(arr));

  let clp = SEED, total = 0;
  for (const [, n, acc] of rows) {
    clp = H([clp, BigInt(acc)]);
    total += parseInt(n, 10);
  }

  console.log("shards        :", rows.length);
  console.log("total variants:", total);
  console.log("C_LP          :", clp.toString());
}

main().then(() => process.exit(0));
