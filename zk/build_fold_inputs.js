// Build the allowed-set Merkle tree and emit fold_inputs.json: a membership path
// for every variant in a key list, plus the running accumulator.
//
//   node build_fold_inputs.js <allowed.bed.gz> <N|all> <outdir> <keys.txt>
//
// keys.txt: one "chrom:pos:ref:alt" per line (1-based pos, chr optional), in the
//           order they will be folded.
//
// Output fold_inputs.json:
//   depth, nLeaves, root, acc_iv, expected_acc,
//   cum_acc[i]  = accumulator after folding variants 0..i,
//   steps[i]    = { chrom, pos, refCode, altCode, pathElements, pathIndices }
//
// acc = 0; acc = Poseidon(acc, Poseidon(chrom,pos,refCode,altCode)) per variant.

const fs = require("fs");
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
  const [bedPath, Narg, outdir, keysPath] = process.argv.slice(2);
  if (!bedPath || !Narg || !outdir || !keysPath) {
    console.error("usage: node build_fold_inputs.js <allowed.bed.gz> <N|all> <outdir> <keys.txt>");
    process.exit(1);
  }

  const keys = fs.readFileSync(keysPath, "utf8")
    .split("\n").map((s) => s.trim().replace(/^chr/i, "")).filter((x) => x.length);
  if (keys.length === 0) throw new Error("empty key list");

  const poseidon = await buildPoseidon();
  const F = poseidon.F;
  const H = (arr) => F.toObject(poseidon(arr));

  const N = Narg === "all" ? Infinity : parseInt(Narg, 10);
  const want = new Set(keys);
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
    if (leaves.length >= N) break;
    if (leaves.length % 1000000 === 0) console.error("  hashed", leaves.length);
  }
  rl.close();
  const nReal = leaves.length;

  const missing = keys.filter((k) => !keyToIdx.has(k));
  if (missing.length) throw new Error(`${missing.length} key(s) not in tree, first: ${missing[0]}`);

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
  const root = levels[levels.length - 1][0];

  const steps = [], cum_acc = [];
  let acc = ACC_IV;
  for (const k of keys) {
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
    cum_acc.push(acc.toString());
    steps.push({
      chrom: chromInt.toString(), pos: posInt.toString(),
      refCode: refCode.toString(), altCode: altCode.toString(),
      pathElements, pathIndices,
    });
  }

  fs.writeFileSync(`${outdir}/fold_inputs.json`, JSON.stringify({
    depth, nLeaves: nReal, root: root.toString(),
    acc_iv: ACC_IV.toString(), expected_acc: cum_acc[cum_acc.length - 1],
    cum_acc, steps,
  }));
  console.error(`depth ${depth}, ${nReal} leaves, ${steps.length} variants; wrote ${outdir}/fold_inputs.json`);
}

main().then(() => process.exit(0));
