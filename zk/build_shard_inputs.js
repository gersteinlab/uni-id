// Parallel variant of build_fold_inputs.js: build the tree once and emit one
// fold-inputs JSON per shard key file, plus each shard's final accumulator.
//
//   node build_shard_inputs.js <allowed.bed.gz> <shards_dir> <out_dir>
//
// shards_dir: shard_0000.keys, shard_0001.keys, ...   (each a key list, in order)
// out_dir   : shard_0000.json, ..., root.txt, and shard_acc.tsv in out_dir/..
//
// Exits nonzero and writes missing_keys.txt if any key is absent from the tree.

const fs = require("fs");
const path = require("path");
const zlib = require("zlib");
const readline = require("readline");
const { buildPoseidon } = require("circomlibjs");

const ACC_IV = 0n;

function chromToInt(c) {
  c = String(c).replace(/^chr/i, "");
  if (c === "X") return 23n;
  if (c === "Y") return 24n;
  if (c === "M" || c === "MT") return 25n;
  return BigInt(c);
}
function alleleToInt(s) {
  const m = { A: 1n, C: 2n, G: 3n, T: 4n, N: 5n };
  let v = 0n;
  for (const ch of String(s).toUpperCase()) {
    const d = m[ch];
    if (d === undefined) throw new Error("bad base '" + ch + "' in '" + s + "'");
    v = v * 6n + d;
  }
  return v;
}

async function main() {
  const [bedPath, shardsDir, outDir] = process.argv.slice(2);
  if (!bedPath || !shardsDir || !outDir) {
    console.error("usage: node build_shard_inputs.js <allowed.bed.gz> <shards_dir> <out_dir>");
    process.exit(1);
  }

  const shardFiles = fs.readdirSync(shardsDir)
    .filter((f) => /^shard_\d+\.keys$/.test(f)).sort();
  if (shardFiles.length === 0) throw new Error("no shard_*.keys in " + shardsDir);

  const shardKeys = {};
  const want = new Set();
  for (const f of shardFiles) {
    const name = f.replace(/\.keys$/, "");
    const ks = fs.readFileSync(path.join(shardsDir, f), "utf8")
      .split("\n").map((s) => s.trim().replace(/^chr/i, "")).filter((x) => x.length);
    shardKeys[name] = ks;
    for (const k of ks) want.add(k);
  }

  const poseidon = await buildPoseidon();
  const F = poseidon.F;
  const H = (arr) => F.toObject(poseidon(arr));

  const keyToIdx = new Map();
  const leaves = [];
  const rl = readline.createInterface({
    input: fs.createReadStream(bedPath).pipe(zlib.createGunzip()),
    crlfDelay: Infinity,
  });
  for await (const line of rl) {
    if (!line || line[0] === "#") continue;
    const f = line.split("\t");
    const chrom = f[0], pos = BigInt(f[1]) + 1n, ref = f[3], alt = f[4];
    const idx = leaves.length;
    leaves.push(H([chromToInt(chrom), pos, alleleToInt(ref), alleleToInt(alt)]));
    const k = `${chrom.replace(/^chr/i, "")}:${pos}:${ref}:${alt}`;
    if (want.has(k)) keyToIdx.set(k, idx);
    if (leaves.length % 1000000 === 0) console.error("  hashed", leaves.length);
  }
  rl.close();
  const nReal = leaves.length;

  const missing = [...want].filter((k) => !keyToIdx.has(k));
  if (missing.length) {
    fs.writeFileSync(path.join(outDir, "missing_keys.txt"), missing.join("\n") + "\n");
    console.error(`${missing.length} key(s) not in tree; see missing_keys.txt`);
    process.exit(2);
  }

  let depth = 0;
  while (1 << depth < nReal) depth++;
  const nPad = 1 << depth;
  const sentinel = H([0n, 0n, 0n, 0n]);
  for (let i = nReal; i < nPad; i++) leaves.push(sentinel);

  const levels = [leaves];
  let cur = leaves;
  while (cur.length > 1) {
    const nx = new Array(cur.length / 2);
    for (let i = 0; i < cur.length; i += 2) nx[i / 2] = H([cur[i], cur[i + 1]]);
    levels.push(nx);
    cur = nx;
  }
  const root = levels[levels.length - 1][0].toString();
  fs.writeFileSync(path.join(outDir, "root.txt"), root + "\n");

  const accLines = [];
  for (const name of Object.keys(shardKeys)) {
    const steps = [];
    let acc = ACC_IV;
    for (const k of shardKeys[name]) {
      const [c, p, r, a] = k.split(":");
      const chromInt = chromToInt(c), posInt = BigInt(p);
      const refCode = alleleToInt(r), altCode = alleleToInt(a);
      let idx = keyToIdx.get(k);
      const pathElements = [], pathIndices = [];
      for (let l = 0; l < depth; l++) {
        const isRight = idx & 1;
        pathElements.push(levels[l][isRight ? idx - 1 : idx + 1].toString());
        pathIndices.push(String(isRight));
        idx = idx >> 1;
      }
      acc = H([acc, H([chromInt, posInt, refCode, altCode])]);
      steps.push({
        chrom: chromInt.toString(), pos: posInt.toString(),
        refCode: refCode.toString(), altCode: altCode.toString(),
        pathElements, pathIndices,
      });
    }
    fs.writeFileSync(path.join(outDir, name + ".json"), JSON.stringify({
      shard: name, depth, nLeaves: nReal, root,
      acc_iv: ACC_IV.toString(), expected_acc: acc.toString(),
      nSteps: steps.length, steps,
    }));
    accLines.push(`${name}\t${steps.length}\t${acc.toString()}`);
  }

  fs.writeFileSync(path.join(outDir, "..", "shard_acc.tsv"), accLines.join("\n") + "\n");
  console.error(`depth ${depth}, ${nReal} leaves, ${shardFiles.length} shards`);
}

main().then(() => process.exit(0));
