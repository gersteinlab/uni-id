# Privacy-stratified partitioning of sequencing data

This toolkit partitions an aligned sequencing file (BAM/CRAM) into a public
layer and two private layers by the allele frequency of the variants each read
carries.

- **pBAM** (public): reference-matched alignments carrying no variant
  information.
- **LP** (less-private): common-variant information.
- **MP** (most-private): rare-variant information, the more sensitive subset.

`pBAM` is the sanitized public file. Adding `LP` recovers common variation;
adding `MP` as well recovers the original file exactly (lossless).

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
- An indexed reference FASTA (`.fai` alongside it)

---

## Pipeline

Three steps: sanitize, createDiff (which also partitions), and restore.
Building the allowed list and error set are one-time inputs, reusable across
samples.

### 1. Sanitize: BAM → pBAM

Produces the public reference-sanitized skeleton.

```bash
sanitize/makepBAM_genome.sh  <input.bam>  <reference.fa>
# writes <input>.sorted.p.bam
```

### 2. createDiff: BAM → LP / MP

Records each read's variant information and routes the read to a layer. A read
goes to LP if every variant it carries is on the allowed list, flagged as a
sequencing error, or a no-call (`N`), or if the read is on a non-primary contig
or below the mapping-quality floor; otherwise it goes to MP. Indels are
left-aligned against the reference before the allowed-list check, so
`--reference` is required for indel-bearing data.

```bash
createDiff/createDiff.py \
    --bam        <input.bam> \
    --region     <chrom> \
    --allowed    allowed_loci.bed.gz \
    --errset     errsites.bed.gz \
    --reference  <reference.fa> \
    --primary-chroms-only \
    --min-mapq   20 \
    --min-bq     0 \
    --gzip \
    --out-prefix <sample>
# writes <sample>.LP.diff.gz and <sample>.MP.diff.gz
```

`createDiff.py` reads a BAM (`--bam`, with `--region` for one contig) or SAM on
stdin (`samtools view <bam> | createDiff.py ...`). Without `--allowed`, every
variant-bearing read routes to MP. Run at base quality 0 (`--min-bq 0`); error
handling is done by the error set, not by base-quality filtering. See
`createDiff.py --help`.

### 3. Restore: pBAM + LP [+ MP] → BAM

Reconstructs a BAM at the chosen access tier from the pBAM and whichever layers
are supplied.

```bash
samtools view -h <sample>.sorted.p.bam \
  | restore/pbam2bam.py --diff <sample>.LP.diff.gz [--diff <sample>.MP.diff.gz] \
                        [--reference <reference.fa>] \
  | samtools view -h -b - > <sample>.restored.bam
```

`--reference` is required only if the supplied layers contain complex
(indel-bearing) reads. Supplying both LP and MP restores the original file
exactly.

---

## Inputs to createDiff

### Allowed list

A bgzipped 5-column BED of common-variant alleles:

```
CHROM   START(0-based)   END(1-based)   REF   ALT
```

The routing keys on `CHROM:END:REF:ALT`. Build one from a directory of
per-chromosome population VCFs (e.g. 1000 Genomes) using the scripts in
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

To exclude trait-associated loci, subtract the positions in the GWAS Catalog
from the allowed list before use; variants at those loci then route to MP. See
`build_gwas_exclusion.sh`.

### Error set

A bgzipped 5-column BED (same format as the allowed list) of substitution sites
that look like sequencing errors: variant allele fraction below a threshold at
adequate depth. Build per chromosome with the `build-errset` mode of
`partition_diff.py`, then concatenate:

```bash
for chr in chr{1..22} chrX; do
  partition/partition_diff.py build-errset \
      --bam <input.bam> --reference <reference.fa> \
      --chrom $chr --vaf-max 0.15 --dp-min 10 \
      --out errsites.$chr.bed.gz
done
zcat errsites.chr*.bed.gz | sort -k1,1 -k2,2n -k4,4 -k5,5 -u | bgzip > errsites.bed.gz
```

The `sort -u` uses the full key (`-k1,1 -k2,2n -k4,4 -k5,5`), not position
alone, so distinct ALT alleles at the same site are not collapsed.

---

## Directory layout

```
sanitize/          BAM -> pBAM              (ptools-derived)
createDiff/        BAM -> LP/MP diffs       (createDiff.py + complex_variants.py)
partition/         build-errset mode        (partition_diff.py)
└── build_allowed/ allowed-list construction
restore/           pBAM + LP[/MP] -> BAM    (pbam2bam.py + complex_variants.py)
```
