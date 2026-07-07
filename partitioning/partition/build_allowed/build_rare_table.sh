#!/bin/bash
# build_rare_table.sh
# Extract variants at or below a MAF threshold from a directory of per-chromosome VCFs.
# Outputs a TSV per chromosome, then concatenates into a genome-wide rare variant table.
#
# Multiallelic sites are split per ALT allele (bcftools norm -m-) and each allele is
# classified on its own frequency. This is the exact complement of build_common_table.sh
# (which keeps MAF > cutoff); together they partition all SNP/indel ALLELES with no
# overlap and no gap at the cutoff. SVs are excluded from both and default to MP.
#
# Usage:
#   bash build_rare_table.sh --vcf-dir <dir> --out-dir <dir> [--maf <float>] [--chr <chr>]
#   bash build_rare_table.sh --out-dir <dir> --merge-only
#
# Arguments:
#   --vcf-dir     Directory containing the 1000G per-chromosome VCF.gz files
#   --out-dir     Output directory for per-chromosome TSVs and final merged table
#   --maf         MAF cutoff (default: 0.05). A variant is RARE if MAF <= cutoff.
#   --chr         Optional: run for a single chromosome only (e.g. chr1).
#   --vcf-pattern Optional: VCF filename pattern; {CHR} is replaced by the chromosome
#                 name. Default matches the 1000G high-coverage panel naming.
#   --merge-only  Skip per-chromosome processing and just merge existing TSVs.
#
# Optional parallelism: if run under a SLURM array (SLURM_ARRAY_TASK_ID set), the task
# id selects the chromosome automatically. This is a convenience only; the script runs
# fine as a plain shell command with --chr, or serially over all chromosomes.
#
# Output:
#   <out-dir>/per_chrom/chr{N}.rare_variants.tsv  - per-chromosome rare variant tables
#   <out-dir>/rare_variants.tsv.gz                - merged, bgzipped genome-wide table
#
# Output TSV columns:
#   CHROM  POS  REF  ALT  AF  MAF  AC  AN  SINGLETON
#   SINGLETON=1 if AC=1 (appears in exactly one haplotype in 1000G)
#
# Requirements: bcftools >= 1.12, bgzip, tabix

set -euo pipefail

# ── Defaults ────────────────────────────────────────────────────────────────
MAF_CUTOFF=0.05
VCF_DIR=""
OUT_DIR=""
SINGLE_CHR=""
MERGE_ONLY=0
# VCF filename pattern; {CHR} is replaced by the chromosome name (e.g. chr1).
VCF_PATTERN="1kGP_high_coverage_Illumina.{CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

# ── Argument parsing ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf-dir)     VCF_DIR="$2"; shift 2 ;;
        --out-dir)     OUT_DIR="$2"; shift 2 ;;
        --maf)         MAF_CUTOFF="$2"; shift 2 ;;
        --chr)         SINGLE_CHR="$2"; shift 2 ;;
        --vcf-pattern) VCF_PATTERN="$2"; shift 2 ;;
        --merge-only)  MERGE_ONLY=1; shift ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ -z "$OUT_DIR" ]]; then
    echo "Error: --out-dir is required."
    exit 1
fi

if [[ "$MERGE_ONLY" -eq 0 && -z "$VCF_DIR" ]]; then
    echo "Error: --vcf-dir is required unless using --merge-only."
    exit 1
fi

mkdir -p "$OUT_DIR/per_chrom"

# ── Merge-only mode ──────────────────────────────────────────────────────────
if [[ "$MERGE_ONLY" -eq 1 ]]; then
    echo "[$(date '+%H:%M:%S')] --merge-only set, merging per-chromosome files from ${OUT_DIR}/per_chrom/..."

    MERGED="${OUT_DIR}/rare_variants.tsv"
    MERGED_GZ="${MERGED}.gz"

    FIRST_TSV=$(ls "${OUT_DIR}/per_chrom/"chr*.rare_variants.tsv 2>/dev/null | head -1)
    if [[ -z "$FIRST_TSV" ]]; then
        echo "Error: no per-chromosome TSV files found in ${OUT_DIR}/per_chrom/"
        exit 1
    fi

    head -1 "$FIRST_TSV" > "$MERGED"
    for CHR in $(seq -f "chr%g" 1 22) chrX; do
        TSV="${OUT_DIR}/per_chrom/${CHR}.rare_variants.tsv"
        if [[ -f "$TSV" ]]; then
            tail -n +2 "$TSV" >> "$MERGED"
        else
            echo "Warning: missing ${TSV}, skipping."
        fi
    done

    echo "[$(date '+%H:%M:%S')] Compressing and indexing..."
    bgzip -f "$MERGED"
    tabix -s 1 -b 2 -e 2 -S 1 "$MERGED_GZ"

    TOTAL=$(zcat "$MERGED_GZ" | tail -n +2 | wc -l)
    SINGLETONS=$(zcat "$MERGED_GZ" | tail -n +2 | awk '$9==1' | wc -l)
    echo ""
    echo "Done. Summary:"
    echo "  Total rare variants (MAF <= ${MAF_CUTOFF}): ${TOTAL}"
    echo "  Of which singletons (AC=1):                 ${SINGLETONS}"
    echo "  Output: ${MERGED_GZ}"
    exit 0
fi

# ── Chromosome list ──────────────────────────────────────────────────────────
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" && -z "$SINGLE_CHR" ]]; then
    TASK_ID="$SLURM_ARRAY_TASK_ID"
    if [[ "$TASK_ID" -le 22 ]]; then
        SINGLE_CHR="chr${TASK_ID}"
    else
        SINGLE_CHR="chrX"
    fi
    echo "SLURM array task ${TASK_ID} → processing ${SINGLE_CHR}"
fi

if [[ -n "$SINGLE_CHR" ]]; then
    CHROMOSOMES=("$SINGLE_CHR")
else
    CHROMOSOMES=($(seq -f "chr%g" 1 22) chrX)
fi

# ── VCF filename pattern ─────────────────────────────────────────────────────
vcf_for_chr() {
    local CHR="$1"
    echo "${VCF_DIR}/${VCF_PATTERN//\{CHR\}/$CHR}"
}

# ── Per-chromosome processing ────────────────────────────────────────────────
for CHR in "${CHROMOSOMES[@]}"; do
    VCF=$(vcf_for_chr "$CHR")
    OUT_TSV="${OUT_DIR}/per_chrom/${CHR}.rare_variants.tsv"

    if [[ ! -f "$VCF" ]]; then
        echo "Warning: VCF not found for ${CHR}: ${VCF} — skipping."
        continue
    fi

    echo "[$(date '+%H:%M:%S')] Processing ${CHR}..."

    # Split multiallelic sites into one record per ALT allele (norm -m-), so each
    # allele is classified independently by its own frequency. AF/AC are Number=A and
    # split correctly. Then restrict to SNP/indel alleles (drops SVs and the spanning-
    # deletion '*' marker). Rare filter keeps MAF = min(AF, 1-AF) <= cutoff; the two-
    # sided form catches rare alleles whether the minor allele is REF or ALT.
    # Note: norm without -f does not left-align indels; add -f <ref.fa> if the panel's
    # indel representations ever need canonicalizing to match the BAM reference.
    bcftools norm -m- "$VCF" \
    | bcftools view \
        --type snps,indels \
        --exclude 'ALT="*"' \
    | bcftools filter \
        --include "INFO/AF <= ${MAF_CUTOFF} || (1 - INFO/AF) <= ${MAF_CUTOFF}" \
    | bcftools query \
        --format '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\t%INFO/AN\n' \
    | awk -v maf_cutoff="$MAF_CUTOFF" '
        BEGIN { OFS="\t"; print "CHROM","POS","REF","ALT","AF","MAF","AC","AN","SINGLETON" }
        {
            chrom=$1; pos=$2; ref=$3; alt=$4; af=$5; ac=$6; an=$7
            maf = (af < 1 - af) ? af : 1 - af
            singleton = (ac == 1) ? 1 : 0
            if (maf <= maf_cutoff)
                print chrom, pos, ref, alt, af, maf, ac, an, singleton
        }
    ' > "$OUT_TSV"

    LINE_COUNT=$(wc -l < "$OUT_TSV")
    echo "[$(date '+%H:%M:%S')] ${CHR}: $((LINE_COUNT - 1)) rare variants written to ${OUT_TSV}"
done

# ── Merge across chromosomes (only when running all at once) ─────────────────
if [[ -z "$SINGLE_CHR" ]]; then
    echo ""
    echo "[$(date '+%H:%M:%S')] Merging per-chromosome files..."

    MERGED="${OUT_DIR}/rare_variants.tsv"
    MERGED_GZ="${MERGED}.gz"

    head -1 "${OUT_DIR}/per_chrom/chr1.rare_variants.tsv" > "$MERGED"
    for CHR in $(seq -f "chr%g" 1 22) chrX; do
        TSV="${OUT_DIR}/per_chrom/${CHR}.rare_variants.tsv"
        if [[ -f "$TSV" ]]; then
            tail -n +2 "$TSV" >> "$MERGED"
        fi
    done

    echo "[$(date '+%H:%M:%S')] Compressing and indexing merged table..."
    bgzip -f "$MERGED"
    tabix -s 1 -b 2 -e 2 -S 1 "$MERGED_GZ"

    TOTAL=$(zcat "$MERGED_GZ" | tail -n +2 | wc -l)
    SINGLETONS=$(zcat "$MERGED_GZ" | tail -n +2 | awk '$9==1' | wc -l)
    echo ""
    echo "Done. Summary:"
    echo "  Total rare variants (MAF <= ${MAF_CUTOFF}): ${TOTAL}"
    echo "  Of which singletons (AC=1):                 ${SINGLETONS}"
    echo "  Output: ${MERGED_GZ}"
fi
