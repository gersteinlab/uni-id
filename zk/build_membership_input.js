// Build the allowed-set Merkle tree and emit input.json for one variant's
// membership proof, plus root.txt.
//
//   node build_membership_input.js <allowed.bed.gz> <N|all> <outdir> [chrom:pos:ref:alt]
//
// allowed.bed.gz: BED with columns chrom, start(0-based), end, ref, alt.
// N|all         : number of leading rows to use, or "all".
// key           : variant to prove (1-based pos, chr prefix optional). If omitted,
//                 the middle leaf is used (requires a finite N).
//
// leaf = Poseidon(chromInt, pos, refCode, altCode); pos = BED start + 1.

const fs = require("fs");
const zlib = require("zlib");
const readline = require("readline");
const { buildPoseidon } = require("circomlibjs");

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
    v = v * 6n + d;                 // base-6, injective over allele strings
  }
  return v;
}

async function main() {
  const bedPath = process.argv[2];
  const Narg = process.argv[3];
  const outdir = process.argv[4];
  let queryKey = process.argv[5] || null;
  if (queryKey) queryKey = queryKey.replace(/^chr/i, "");
  if (!bedPath || !Narg || !outdir) {
    console.error("usage: node build_membership_input.js <allowed.bed.gz> <N|all> <outdir> [chrom:pos:ref:alt]");
    process.exit(1);
  }

  const poseidon = await buildPoseidon();
  const F = poseidon.F;
  const H = (arr) => F.toObject(poseidon(arr));

  const N = Narg === "all" ? Infinity : parseInt(Narg, 10);
  const middleIdx = Number.isFinite(N) ? Math.floor(N / 2) : null;
  const doMap = !!queryKey;
  const keyToIdx = doMap ? new Map() : null;
  let qFields = null;

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
    if (doMap) keyToIdx.set(`${chrom.replace(/^chr/i, "")}:${pos}:${ref}:${alt}`, idx);
    if (idx === middleIdx) qFields = { chrom, pos: pos.toString(), ref, alt };
    if (leaves.length >= N) break;
    if (leaves.length % 1000000 === 0) console.error("  hashed", leaves.length);
  }
  rl.close();
  const nReal = leaves.length;
  if (nReal === 0) throw new Error("no leaves read");

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

  let qIdx;
  if (queryKey) {
    if (!keyToIdx.has(queryKey)) throw new Error("query key not in tree: " + queryKey);
    qIdx = keyToIdx.get(queryKey);
    const [c, p, r, a] = queryKey.split(":");
    qFields = { chrom: c, pos: p, ref: r, alt: a };
  } else {
    if (middleIdx === null) throw new Error("no query key given and N is 'all'");
    qIdx = middleIdx;
  }

  const pathElements = [], pathIndices = [];
  let idx = qIdx;
  for (let l = 0; l < depth; l++) {
    const isRight = idx & 1;
    pathElements.push(levels[l][isRight ? idx - 1 : idx + 1].toString());
    pathIndices.push(String(isRight));
    idx = idx >> 1;
  }

  const input = {
    chrom: chromToInt(qFields.chrom).toString(),
    pos: qFields.pos,
    refCode: alleleToInt(qFields.ref).toString(),
    altCode: alleleToInt(qFields.alt).toString(),
    pathElements,
    pathIndices,
    root: root.toString(),
  };
  fs.writeFileSync(outdir + "/input.json", JSON.stringify(input, null, 2));
  fs.writeFileSync(outdir + "/root.txt", root.toString() + "\n");
  console.error(`depth ${depth}, ${nReal} leaves; wrote ${outdir}/input.json and root.txt`);
}

main().then(() => process.exit(0));
