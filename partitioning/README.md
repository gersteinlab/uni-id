# uni-id: privacy-stratified partitioning of sequencing data

This toolkit partitions an aligned sequencing file (BAM/CRAM) into privacy
layers, so that broadly shareable content can be released while
identifying content is held back behind stricter access control.

Four layers result:

- **pBAM** (public): a reference-sanitized skeleton in which every read's
  sequence is replaced by the reference. Shareable without disclosing variants.
- **LP** (less-private): the edits needed to restore reads whose variation is
  *non-characterizing* — common variants, sequencing errors, no-calls, and
  reads in unreliable regions.
- **MP** (most-private): the edits needed to restore reads carrying rare or
  novel (identifying) variation.

Restoring `pBAM` alone gives the public skeleton; `pBAM + LP` restores common
variation; `pBAM + LP + MP` restores the original file exactly (lossless).

This work extends [ptools](https://github.com/) (Gürsoy et al., *Cell* 2020);
files derived from ptools carry an attribution header and remain under the MIT
License. See `LICENSE`.

---

## Dependencies

- Python 3 (standard library only for the core scripts)
- [samtools](http://www.htslib.org/) (≥ 1.11)
- [bcftools](http://www.htslib.org/) (≥ 1.12) — for building the error set and allowed list
- [Picard](https://broadinstitute.github.io/picard/) 2.23.8 — for the sanitize step
- [plink](https://www.cog-genomics.org/plink/) 1.9 — for building the allowed list
- `bgzip` / `tabix` (htslib)

---

## The four operations

The pipeline is four steps. Each is a single command; only the allowed-list and
error-set construction (inputs to partitioning) take more than one line, and
those are one-time or reusable across many samples.

### 1. Sanitize: BAM → pBAM

Produces the public reference-sanitized skeleton.

```bash
sanitize/makepBAM_genome.sh  <input.bam>  <reference.fa>
# writes <input>.sorted.p.bam
```

### 2. createDiff: BAM → diff

Records, for every read, the edits needed to restore it from the pBAM. Run at
base quality 0 (`--min-bq 0`): error handling is done later by the partition
step via allele fraction, not by base-quality filtering here.

```bash
createDiff/createDiff.py --bam <input.bam> --min-bq 0 > <sample>.diff
```

`createDiff.py` reads a BAM directly (`--bam`, enabling QNAME-drop) or SAM on
stdin (`samtools view <bam> | createDiff.py`). See `createDiff.py --help`.

### 3. Partition: diff → LP / MP

Splits the diff into the two private layers. A read goes to LP only if *every*
substitution it carries is non-characterizing: on the allowed list, a
sequencing error, a no-call (`N`), inside a blacklist region, or on a read that
is on a non-primary contig or below a MAPQ threshold. Everything else goes to MP.

Partitioning needs two inputs described below — an **allowed list** (common
variants) and, optionally, an **error set** (low-VAF sites). It can also take a
blacklist and MAPQ / primary-contig rules.

```bash
partition/partition_diff.py partition \
    --diff       <sample>.diff \
    --allowed    allowed_loci.bed.gz \
    --errset     errsites.bed.gz \
    --blacklist  hg38-blacklist.v2.bed.gz \
    --primary-chroms-only \
    --min-mapq   30 \
    --gzip \
    --out-prefix <sample>
# writes <sample>.LP.diff.gz and <sample>.MP.diff.gz
```

Only `--diff`, `--allowed`, and `--out-prefix` are required; the rest are
optional filters. See `partition_diff.py partition --help`.

### 4. Restore: pBAM + LP [+ MP] → BAM

Reconstructs a BAM at the chosen access tier from the pBAM and whichever layers
are supplied.

```bash
samtools view -h <sample>.sorted.p.bam \
  | restore/pbam2bam.py --diff <sample>.LP.diff.gz [--diff <sample>.MP.diff.gz] \
                        [--reference <reference.fa>] \
  | samtools view -h -b - > <sample>.restored.bam
```

`--reference` is required only if the supplied layers contain structured
complex (indel-bearing) reads. Supplying both LP and MP restores the original
file exactly.

---

## Inputs to the partition step

### Allowed list

A bgzipped 5-column BED of common-variant alleles:

```
CHROM   START(0-based)   END(1-based)   REF   ALT
```

The partition step keys on `CHROM:END:REF:ALT`. Build one from a directory of
per-chromosome population VCFs (e.g. 1000 Genomes) using the three scripts in
`partition/build_allowed/`:

```bash
# 1. common variants (MAF > cutoff) and rare variants (MAF <= cutoff)
partition/build_allowed/build_common_table.sh --vcf-dir <vcfs> --out-dir <out> --maf 0.05
partition/build_allowed/build_rare_table.sh   --vcf-dir <vcfs> --out-dir <out> --maf 0.05

# 2. allowed list = common alleles NOT in LD with any rare allele
partition/build_allowed/build_allowed_list.sh \
    --vcf-dir      <vcfs> \
    --rare-table   <out>/rare_variants.tsv.gz \
    --common-table <out>/common_variants.tsv.gz \
    --out-dir      <out>
# writes <out>/allowed_loci.bed.gz
```

Each script takes `--vcf-pattern` if your VCF filenames differ from the 1000G
default, and `--chr` to run one chromosome at a time (parallelizable; a SLURM
array is detected automatically if present, but is optional).

### Error set (optional)

A bgzipped 5-column BED (same format as the allowed list) of substitution sites
that look like sequencing errors: variant allele fraction below a threshold at
adequate depth. These are treated as non-characterizing (routed to LP). Build
per chromosome with the `build-errset` mode of `partition_diff.py`, then
concatenate:

```bash
for chr in chr{1..22} chrX; do
  partition/partition_diff.py build-errset \
      --bam <input.bam> --reference <reference.fa> \
      --chrom $chr --vaf-max 0.15 --dp-min 10 \
      --out errsites.$chr.bed.gz
done
zcat errsites.chr*.bed.gz | sort -k1,1 -k2,2n -k4,4 -k5,5 -u | bgzip > errsites.bed.gz
```

Note the `sort -u` uses the full key (`-k1,1 -k2,2n -k4,4 -k5,5`), not just
position — deduplicating on position alone would silently drop distinct ALT
alleles at the same site.

### Blacklist (optional)

Reads whose substitutions fall in problem regions (repetitive, low-mappability)
are treated as non-characterizing. This toolkit does not ship a blacklist; the
ENCODE hg38 unified blacklist is available from the Boyle Lab:

<https://github.com/Boyle-Lab/Blacklist> (Amemiya et al., *Sci Rep* 2019).

---

## Directory layout

```
partitioning/
├── sanitize/          BAM -> pBAM              (ptools-derived)
├── createDiff/        BAM -> diff              (createDiff.py + complex_variants.py)
├── partition/         diff -> LP/MP            (partition_diff.py; build-errset mode)
│   └── build_allowed/ allowed-list construction (3 scripts)
└── restore/           pBAM + LP[/MP] -> BAM    (pbam2bam.py + complex_variants.py)
```
