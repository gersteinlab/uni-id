#!/bin/bash
# build_allowed_list.sh
# Construct the partition "allowed list": common variant alleles that are NOT in LD
# with any rare variant. A read whose every carried allele is on this list goes to LP;
# anything else (rare, LD-tainted common, multiallelic alt not listed, novel, SV) goes
# to MP by the fail-closed default.
#
#   allowed = (common alleles) - (common alleles in LD r2 >= R2 within WINDOW of a rare allele)
#
# Per chromosome:
#   1. Build a normalized, allele-split, ID-annotated VCF (chr:pos:ref:alt IDs) and a
#      plink binary from it. These are MAF-INDEPENDENT and cached, so a multi-cutoff
#      sweep reuses them — only the rare/common lists below change with the cutoff.
#   2. Take the rare table as plink --ld-snp-list anchors.
#   3. plink --r2 reports common partners in LD with a rare anchor -> "tainted" commons.
#   4. allowed_chrom = (common alleles on chrom) - (tainted commons).
#   5. Emit a 5-col BED: CHROM  START(0-based)  END  REF  ALT.
# Merge concatenates, sorts, bgzips, tabix-indexes -> allowed_loci.bed.gz
#
# Matching is by the bcftools-set ID string (chr:pos:ref:alt), so plink's habit of
# rewriting the CHROM column to bare numbers does not affect anything here.
#
# Usage:
#   bash build_allowed_list.sh --vcf-dir <dir> --rare-table <rare.tsv.gz> \
#        --common-table <common.tsv.gz> --out-dir <dir> \
#        [--cache-dir <dir>] [--work-dir <dir>] \
#        [--r2 <float>] [--window-kb <int>] [--chr <chr>]
#   bash build_allowed_list.sh --out-dir <dir> --merge-only
#
# Arguments:
#   --vcf-dir       Directory with the per-chromosome VCF.gz files
#   --rare-table    Merged rare_variants.tsv.gz (LD anchors; defines the cutoff)
#   --common-table  Merged common_variants.tsv.gz (allowed-list candidates)
#   --out-dir       Output dir for per-chrom BEDs and the final merged BED
#   --cache-dir     MAF-independent plink binaries (default: <out-dir>/cache).
#                   Share this across MAF cutoffs to avoid rebuilding genotypes.
#   --work-dir      Scratch for large intermediate .ld files (default: <out-dir>/work).
#                   These can be tens of GB per chromosome and are deleted after use.
#   --r2            LD r^2 threshold (default: 0.25)
#   --window-kb     LD window in kb (default: 500)
#   --vcf-pattern   VCF filename pattern; {CHR} replaced by chromosome (default: 1000G naming)
#   --chr           Single chromosome (e.g. chr1); optionally set automatically under SLURM array
#   --merge-only    Just merge existing per-chrom BEDs
#
# Requirements: bcftools >= 1.12, plink 1.9, bgzip, tabix

set -euo pipefail

# ── Defaults ────────────────────────────────────────────────────────────────
# VCF filename pattern; {CHR} is replaced by the chromosome name (e.g. chr1).
VCF_PATTERN="1kGP_high_coverage_Illumina.{CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
VCF_DIR=""
RARE_TABLE=""
COMMON_TABLE=""
OUT_DIR=""
CACHE_DIR=""
WORK_DIR=""
R2=0.25
WINDOW_KB=500
SINGLE_CHR=""
MERGE_ONLY=0
MEM_MB=40000        # plink workspace cap (MB); must sit under the SBATCH --mem cgroup
THREADS=8           # plink threads; match SBATCH --cpus-per-task
KEEP=""             # optional plink --keep file (e.g. unrelated_2504.keep); applied to --r2 only
CACHE_ONLY=0        # if 1, build only the normalized VCF + plink binary cache, then exit

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf-dir)       VCF_DIR="$2"; shift 2 ;;
        --rare-table)    RARE_TABLE="$2"; shift 2 ;;
        --common-table)  COMMON_TABLE="$2"; shift 2 ;;
        --out-dir)       OUT_DIR="$2"; shift 2 ;;
        --cache-dir)     CACHE_DIR="$2"; shift 2 ;;
        --work-dir)      WORK_DIR="$2"; shift 2 ;;
        --r2)            R2="$2"; shift 2 ;;
        --window-kb)     WINDOW_KB="$2"; shift 2 ;;
        --vcf-pattern)   VCF_PATTERN="$2"; shift 2 ;;
        --mem-mb)        MEM_MB="$2"; shift 2 ;;
        --threads)       THREADS="$2"; shift 2 ;;
        --keep)          KEEP="$2"; shift 2 ;;
        --cache-only)    CACHE_ONLY=1; shift ;;
        --chr)           SINGLE_CHR="$2"; shift 2 ;;
        --merge-only)    MERGE_ONLY=1; shift ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ -z "$OUT_DIR" ]]; then
    echo "Error: --out-dir is required."
    exit 1
fi
CACHE_DIR="${CACHE_DIR:-$OUT_DIR/cache}"
WORK_DIR="${WORK_DIR:-$OUT_DIR/work}"
mkdir -p "$OUT_DIR/per_chrom" "$CACHE_DIR" "$WORK_DIR"

# ── Merge-only mode ──────────────────────────────────────────────────────────
if [[ "$MERGE_ONLY" -eq 1 ]]; then
    echo "[$(date '+%H:%M:%S')] --merge-only: merging per-chrom BEDs from ${OUT_DIR}/per_chrom/..."
    MERGED_GZ="${OUT_DIR}/allowed_loci.bed.gz"

    if ! ls "${OUT_DIR}/per_chrom/"chr*.allowed.bed >/dev/null 2>&1; then
        echo "Error: no per-chrom BED files found in ${OUT_DIR}/per_chrom/"
        exit 1
    fi

    cat "${OUT_DIR}/per_chrom/"chr*.allowed.bed \
        | sort -k1,1 -k2,2n \
        | bgzip > "$MERGED_GZ"
    tabix -p bed "$MERGED_GZ"

    TOTAL=$(zcat "$MERGED_GZ" | wc -l)
    echo ""
    echo "Done. Summary:"
    echo "  Total allowed (common, non-LD-tainted) alleles: ${TOTAL}"
    echo "  Output: ${MERGED_GZ}"
    exit 0
fi

# ── Required args for per-chrom mode ─────────────────────────────────────────
# In cache-only mode only the VCF is needed (rare/common tables are not used).
if [[ "$CACHE_ONLY" -eq 1 ]]; then
    REQUIRED=(VCF_DIR)
else
    REQUIRED=(VCF_DIR RARE_TABLE COMMON_TABLE)
fi
for req in "${REQUIRED[@]}"; do
    if [[ -z "${!req}" ]]; then
        echo "Error: --${req,,} is required (unless --merge-only)." | tr '_' '-'
        exit 1
    fi
done

# ── Chromosome selection ─────────────────────────────────────────────────────
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" && -z "$SINGLE_CHR" ]]; then
    TASK_ID="$SLURM_ARRAY_TASK_ID"
    if [[ "$TASK_ID" -le 22 ]]; then SINGLE_CHR="chr${TASK_ID}"; else SINGLE_CHR="chrX"; fi
    echo "SLURM array task ${TASK_ID} → processing ${SINGLE_CHR}"
fi
if [[ -z "$SINGLE_CHR" ]]; then
    echo "Error: per-chrom mode needs --chr or a SLURM array context."
    exit 1
fi
CHR="$SINGLE_CHR"

VCF="${VCF_DIR}/${VCF_PATTERN//\{CHR\}/$CHR}"
if [[ ! -f "$VCF" ]]; then
    echo "Error: VCF not found for ${CHR}: ${VCF}"
    exit 1
fi

NORMVCF="${CACHE_DIR}/${CHR}.ids.vcf.gz"
BEDPREFIX="${CACHE_DIR}/${CHR}"
RARE_SNPS="${WORK_DIR}/${CHR}.rare.snps"
COMMON_IDS="${WORK_DIR}/${CHR}.common.ids"
LD_PREFIX="${WORK_DIR}/${CHR}.ld"
TAINTED="${WORK_DIR}/${CHR}.tainted.ids"
ALLOWED_IDS="${WORK_DIR}/${CHR}.allowed.ids"
OUT_BED="${OUT_DIR}/per_chrom/${CHR}.allowed.bed"

# ── 1. Normalized, allele-split, ID-annotated VCF (cached, MAF-independent) ──
if [[ ! -f "$NORMVCF" ]]; then
    echo "[$(date '+%H:%M:%S')] ${CHR}: building normalized ID-annotated VCF..."
    bcftools norm -m- "$VCF" \
        | bcftools view --type snps,indels --exclude 'ALT="*"' \
        | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o "$NORMVCF"
    tabix -p vcf "$NORMVCF"
else
    echo "[$(date '+%H:%M:%S')] ${CHR}: reusing cached ${NORMVCF}"
fi

# ── 2. plink binary (cached, MAF-independent) ───────────────────────────────
# IDs come from bcftools, so --keep-allele-order just preserves REF/ALT; matching is
# ID-based regardless of how plink rewrites the CHROM column. --allow-extra-chr tolerates
# the 'chr' prefix / chrX without crashing.
if [[ ! -f "${BEDPREFIX}.bed" ]]; then
    echo "[$(date '+%H:%M:%S')] ${CHR}: building plink binary..."
    plink --vcf "$NORMVCF" \
          --double-id \
          --keep-allele-order \
          --allow-extra-chr \
          --memory "$MEM_MB" \
          --threads "$THREADS" \
          --make-bed \
          --out "$BEDPREFIX"
else
    echo "[$(date '+%H:%M:%S')] ${CHR}: reusing cached ${BEDPREFIX}.bed"
fi

if [[ "$CACHE_ONLY" -eq 1 ]]; then
    echo "[$(date '+%H:%M:%S')] ${CHR}: cache-only mode, done (norm VCF + plink binary ready)."
    exit 0
fi

# ── 3. Rare anchor list (MAF-dependent) ─────────────────────────────────────
echo "[$(date '+%H:%M:%S')] ${CHR}: extracting rare anchors..."
zcat "$RARE_TABLE" | awk -v c="$CHR" 'NR>1 && $1==c {print $1":"$2":"$3":"$4}' > "$RARE_SNPS"
N_RARE=$(wc -l < "$RARE_SNPS")
echo "[$(date '+%H:%M:%S')] ${CHR}: ${N_RARE} rare anchors."

# ── 4. Common candidate IDs (MAF-dependent), sorted unique ──────────────────
zcat "$COMMON_TABLE" | awk -v c="$CHR" 'NR>1 && $1==c {print $1":"$2":"$3":"$4}' \
    | sort -u > "$COMMON_IDS"
N_COMMON=$(wc -l < "$COMMON_IDS")
echo "[$(date '+%H:%M:%S')] ${CHR}: ${N_COMMON} common candidates."

# ── 5. plink LD: report partners of rare anchors at r2 >= R2 within WINDOW ───
# --keep (if given) restricts LD computation to a sample subset, e.g. the 2,504
# unrelated samples, so r^2 reflects population LD rather than trio relatedness.
# Applied here only — the cached binary keeps all 3,202 so the subset stays a knob.
KEEP_FLAG=()
if [[ -n "$KEEP" ]]; then
    if [[ ! -f "$KEEP" ]]; then
        echo "Error: --keep file not found: ${KEEP}"
        exit 1
    fi
    KEEP_FLAG=(--keep "$KEEP")
    echo "[$(date '+%H:%M:%S')] ${CHR}: restricting LD to samples in ${KEEP}"
fi

echo "[$(date '+%H:%M:%S')] ${CHR}: running plink --r2 (r2>=${R2}, ${WINDOW_KB}kb)..."
plink --bfile "$BEDPREFIX" \
      --allow-extra-chr \
      --memory "$MEM_MB" \
      --threads "$THREADS" \
      "${KEEP_FLAG[@]}" \
      --r2 \
      --ld-snp-list "$RARE_SNPS" \
      --ld-window-kb "$WINDOW_KB" \
      --ld-window 99999 \
      --ld-window-r2 "$R2" \
      --out "$LD_PREFIX"

# ── 6. Tainted commons = common IDs appearing as an LD partner (SNP_B, col 6) ─
# Stream the (potentially large) .ld once, keep only partners that are common, then
# delete the .ld immediately to reclaim disk.
echo "[$(date '+%H:%M:%S')] ${CHR}: extracting LD-tainted common alleles..."
awk 'NR==FNR { c[$1]=1; next } FNR>1 && ($6 in c) { print $6 }' \
    "$COMMON_IDS" "${LD_PREFIX}.ld" \
    | sort -u > "$TAINTED"
N_TAINTED=$(wc -l < "$TAINTED")
rm -f "${LD_PREFIX}.ld"
echo "[$(date '+%H:%M:%S')] ${CHR}: ${N_TAINTED} common alleles tainted by LD with a rare allele."

# ── 7. allowed = common - tainted, emit BED ─────────────────────────────────
comm -23 "$COMMON_IDS" "$TAINTED" > "$ALLOWED_IDS"
N_ALLOWED=$(wc -l < "$ALLOWED_IDS")

# ID is chr:pos:ref:alt — chrom and pos contain no ':', so split is unambiguous.
awk -F: 'BEGIN{OFS="\t"} { print $1, $2-1, $2, $3, $4 }' "$ALLOWED_IDS" \
    | sort -k1,1 -k2,2n > "$OUT_BED"

echo "[$(date '+%H:%M:%S')] ${CHR}: ${N_ALLOWED} allowed alleles → ${OUT_BED}"
echo "[$(date '+%H:%M:%S')] ${CHR}: summary  common=${N_COMMON}  tainted=${N_TAINTED}  allowed=${N_ALLOWED}"
