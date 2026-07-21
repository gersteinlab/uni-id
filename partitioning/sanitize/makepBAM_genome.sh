#!/bin/bash
# ----------------------------------------------------------------------------
# This file is part of the uni-id partitioning toolkit.
#
# Portions derived from ptools (Gürsoy et al., "Data Sanitization to Reduce
# Private Information Leakage from Functional Genomics", Cell 2020), which is
# distributed under the MIT License. The original ptools copyright and MIT
# permission notice are retained below; modifications for uni-id are also
# released under the MIT License.
#
# MIT License. Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation files, to deal
# in the Software without restriction, including the rights to use, copy,
# modify, merge, publish, distribute, sublicense, and/or sell copies. The above
# notice and this permission notice shall be included in all copies. THE
# SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.
# ----------------------------------------------------------------------------
#
# makepBAM_genome.sh -- sanitize a BAM into a reference-matched pBAM.
#
# Reads are split on whether their CIGAR contains N (spliced alignments need
# different handling), each half is rewritten to match the reference, and the
# results are merged and sorted.
#
# Usage:
#   makepBAM_genome.sh <input.bam> <reference.fa> [options]
#
# Options:
#   --region R     restrict to region R. Repeatable. Requires an indexed BAM.
#                  Contig names containing colons (HLA alleles, e.g.
#                  HLA-A*01:01:01:01) MUST use samtools brace syntax:
#                      --region '{HLA-A*01:01:01:01}'
#                  Without braces, samtools reads the name as contig:start-end
#                  and silently returns the wrong reads.
#   --out PREFIX   output prefix (default: input basename minus .bam)
#   --tmp DIR      scratch directory (default: beside the output)
#   --threads N    samtools threads (default 4)
#
# Writes <PREFIX>.sorted.p.bam
#
# Requires: samtools, numpy (getSeq_genome_*.py), Biopython (PrintSequence.py).
# Unlike createDiff/partition/restore, this step is NOT standard-library-only.
#
# Scale note: the intermediate SAM is uncompressed and roughly 5x the BAM. A
# whole 30x genome is ~200 GB of scratch and runs single-threaded through
# python. For genome-scale input, parallelize with --region (one job per
# chromosome) and `samtools merge` the results. Whatever scheme you use, make
# sure it covers EVERY contig in the BAM header: a pBAM missing contigs cannot
# reconstruct the diff records that reference them.
# ----------------------------------------------------------------------------
set -eu

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

usage() { sed -n '18,45p' "${BASH_SOURCE[0]}" | sed 's/^# \?//'; exit 1; }

[[ $# -lt 2 ]] && usage
BAM=$1; REF=$2; shift 2

REGIONS=()
OUTPREFIX=""
TMPDIR_ARG=""
THREADS=4
while [[ $# -gt 0 ]]; do
  case "$1" in
    --region)  REGIONS+=("$2"); shift 2 ;;
    --out)     OUTPREFIX="$2"; shift 2 ;;
    --tmp)     TMPDIR_ARG="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "unknown option: $1"; usage ;;
  esac
done

command -v samtools >/dev/null || { echo "samtools not on PATH"; exit 1; }
# getSeq_genome_*.py need numpy, and PrintSequence.py needs Biopython. Check
# BEFORE doing any work: the split and extraction can take an hour, and
# discovering a missing import after that wastes the whole run.
python3 -c "import numpy" 2>/dev/null || {
  echo "ERROR: numpy not available to python3 (needed by getSeq_genome_*.py)"
  echo "  activate an environment that has it, e.g.:"
  echo "    module load miniconda && conda activate ptools"
  exit 1; }
python3 -c "from Bio import SeqIO" 2>/dev/null || {
  echo "ERROR: Biopython not available to python3 (needed by PrintSequence.py)"
  echo "  activate an environment that has it, e.g.:"
  echo "    module load miniconda && conda activate ptools"
  exit 1; }
[[ -s "$BAM" ]] || { echo "missing input: $BAM"; exit 1; }
[[ -s "$REF" ]] || { echo "missing reference: $REF"; exit 1; }
# resolved relative to this script, not $PATH: a published package should not
# require its own helpers to be installed somewhere global
for p in getSeq_genome_wN.py getSeq_genome_woN.py; do
  [[ -s "$HERE/$p" ]] || { echo "missing $HERE/$p"; exit 1; }
done

if [[ -z "$OUTPREFIX" ]]; then
  b=$(basename "$BAM"); OUTPREFIX=${b%.bam}
fi
TMP=${TMPDIR_ARG:-$(dirname "$OUTPREFIX")/pbam_tmp_$$}
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

if [[ ${#REGIONS[@]} -gt 0 ]]; then
  [[ -s "${BAM}.bai" || -s "${BAM}.csi" ]] || {
    echo "--region needs an indexed BAM; run: samtools index $BAM"; exit 1; }
  echo "regions: ${REGIONS[*]}"
fi

# ---- one pass; split unplaced off, then split the rest on N in the CIGAR -----
# The original read the BAM twice, once per half. One pass, three streams.
#
# Unplaced reads (RNAME '*') must not reach getSeq_genome_*.py: those look the
# reference up at RNAME/POS, and an unplaced read has neither. There is no
# reference at a nonexistent position, so the only sanitization available is
# masking SEQ to N. Their bases are exactly as identifying as any other read's,
# so leaving them in the clear would be a leak.
echo "[$(date '+%H:%M:%S')] splitting on N-CIGAR, separating unplaced reads"
samtools view -H "$BAM" > "$TMP/header.txt"
samtools view -@ "$THREADS" "$BAM" ${REGIONS[@]+"${REGIONS[@]}"} \
  | awk -v a="$TMP/withN.sam" -v b="$TMP/withoutN.sam" -v u="$TMP/unplaced.sam" \
        '{ if ($3 == "*")      print > u
           else if ($6 ~ /N/)  print > a
           else                print > b }'
touch "$TMP/withN.sam" "$TMP/withoutN.sam" "$TMP/unplaced.sam"  # awk writes nothing for an empty stream
echo "  with N:    $(wc -l < "$TMP/withN.sam")"
echo "  without N: $(wc -l < "$TMP/withoutN.sam")"
echo "  unplaced:  $(wc -l < "$TMP/unplaced.sam")"

if [[ ! -s "$TMP/withN.sam" && ! -s "$TMP/withoutN.sam" && ! -s "$TMP/unplaced.sam" ]]; then
  echo "no reads selected; nothing to do"; exit 1
fi

# ---- sanitize each half ------------------------------------------------------
PARTS=()
if [[ -s "$TMP/withN.sam" ]]; then
  echo "[$(date '+%H:%M:%S')] sanitizing spliced reads"
  python3 "$HERE/getSeq_genome_wN.py" "$REF" "$TMP/header.txt" "$TMP/withN.sam" \
    | samtools view -h -b -@ "$THREADS" - > "$TMP/withN.p.bam"
  PARTS+=("$TMP/withN.p.bam")
fi
if [[ -s "$TMP/withoutN.sam" ]]; then
  echo "[$(date '+%H:%M:%S')] sanitizing unspliced reads"
  python3 "$HERE/getSeq_genome_woN.py" "$REF" "$TMP/header.txt" "$TMP/withoutN.sam" \
    | samtools view -h -b -@ "$THREADS" - > "$TMP/withoutN.p.bam"
  PARTS+=("$TMP/withoutN.p.bam")
fi
if [[ -s "$TMP/unplaced.sam" ]]; then
  echo "[$(date '+%H:%M:%S')] sanitizing unplaced reads (masking SEQ to N)"
  # QUAL is kept, as in every other pBAM: it is invariant under sanitization
  # and the diff recovers it from here rather than storing it.
  { cat "$TMP/header.txt"
    awk -F'\t' 'BEGIN{OFS="\t"}
         { n = length($10); s = ""
           for (i = 0; i < n; i++) s = s "N"
           $10 = s
           print }' "$TMP/unplaced.sam"
  } | samtools view -h -b -@ "$THREADS" - > "$TMP/unplaced.p.bam"
  PARTS+=("$TMP/unplaced.p.bam")
fi

# ---- merge + sort ------------------------------------------------------------
echo "[$(date '+%H:%M:%S')] merging and sorting"
if [[ ${#PARTS[@]} -eq 1 ]]; then
  cp "${PARTS[0]}" "$TMP/merged.p.bam"
else
  samtools merge -f -@ "$THREADS" "$TMP/merged.p.bam" "${PARTS[@]}"
fi
# samtools sort puts unplaced reads (RNAME '*') last, which is where a
# coordinate-sorted BAM already keeps them, so pBAM line N stays BAM line N.
# That correspondence is what lets the diff index records by line number
# instead of a content key.
samtools sort -@ "$THREADS" -T "$TMP/sort" "$TMP/merged.p.bam" \
  -o "${OUTPREFIX}.sorted.p.bam"
samtools index -@ "$THREADS" "${OUTPREFIX}.sorted.p.bam"

# a truncated BAM indexes without complaint and fails later, quietly
samtools quickcheck -v "${OUTPREFIX}.sorted.p.bam" || { echo "output BAM CORRUPT"; exit 1; }
echo "[$(date '+%H:%M:%S')] done: ${OUTPREFIX}.sorted.p.bam"
echo "  reads: $(samtools view -c -@ "$THREADS" "${OUTPREFIX}.sorted.p.bam")"
echo "  size:  $(du -h "${OUTPREFIX}.sorted.p.bam" | cut -f1)"
