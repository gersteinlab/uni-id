#!/usr/bin/env python3
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
# createDiff.py
#
# Reads an aligned BAM and writes the two private layers directly:
#
#     BAM  ->  LP.diff  +  MP.diff
#
# Each record names its read by POSITION IN THE pBAM, not by a content key.
# That is sound because pBAM line N is BAM line N: the sanitize transform
# preserves QNAME, FLAG, RNAME, POS, MAPQ, RNEXT, PNEXT, TLEN and QUAL, and
# alters only CIGAR. Measured on NA12878: chr20 (17,705,654 reads) and chr1
# (63,275,059 reads), 100.0000% agreement on both read identity and POS.
#
# Why not a content key. The previous format identified a read by columns
# 0-8. Measured on the chr6 pBAM, that key has 52 duplicates, because
# supplementary alignments of one chimeric read share QNAME, FLAG and all mate
# fields by construction. Two reads sharing a key means one record serves both
# and the other is restored from the wrong record. No subset of invariant
# fields separates them, so the key cannot be repaired. A line number can
# never collide.
#
# Why classification happens here. Routing needs RNAME, MAPQ and the read's
# variants. Reading them from the BAM as it streams means they never have to be
# stored, which removes RNAME, MAPQ and the whole 9-column key from the format,
# and removes the intermediate whole-genome diff entirely.
#
# Output format (tab-separated). RESTORE is the trailing optional-tag columns.
#
#   S:  DELTA  S  subs                                RESTORE...
#   C:  DELTA  C  CIGAR  variants  splices  forcemp   RESTORE...
#   U:  DELTA  U  seq                                 RESTORE...
#
#   DELTA     pBAM lines advanced since the previous record IN THIS FILE. The
#             first record carries the absolute line number. Delta rather than
#             absolute because absolute numbers are ~9 digits and all differ,
#             so they barely compress; gaps average ~4 and compress to ~1 byte.
#   subs      "gpos:ref>alt;..." genomic coordinates, so the allowed list can
#             be tested with no offset arithmetic.
#   CIGAR     stored ONLY for C: it is the one field sanitize destroys.
#   variants  "gpos:ref>alt;..." substitutions AND indels.
#   splices   "off:I:bases;off:S:bases;..." inserted / soft-clipped bases.
#   forcemp   "1" if the read must stay MP (soft clip / unanchorable indel).
#   seq       raw SEQ, for reads that cannot be rebuilt from the reference.
#
# QNAME, FLAG, POS, RNEXT, PNEXT, TLEN and QUAL are not stored: all are
# recovered from the pBAM at line N.
#
# Reference-only reads get no record at all, and there is no flag to change
# that. The old --keep-ref-only put them in LP, which makes not-LP EXACTLY
# equal to MP and hands an LP holder a perfect list of every read carrying a
# rare variant. Dropping them is what makes the layering private.
#
# They are the ~71% of reads that
# match the reference (548,094,823 of 768,580,569 in NA12878), and they
# reconstruct by passing the pBAM line through untouched. This is also what
# makes the layering private: an LP holder can tell which reads are NOT LP, but
# that set is MP plus ref-only, and MP is only ~7.6% of it.
#
# Granularity is a property of invocation, not of the format. Whole BAM in ->
# one pBAM, one LP/MP pair, global numbering. --region chr7 in -> the chr7
# pBAM, a chr7 LP/MP pair, chr7-local numbering. Sanitize the same region with
# the same flag and the pairing holds. Never mix a per-region diff with a
# whole-genome pBAM.
#
# --allowed is optional. Without it the membership set is empty, so no
# difference is non-characterizing and every variant-bearing read goes to MP.
# That is the fail-closed reading of the layers rather than a degenerate case,
# and it avoids loading ~13 GB of allowed list for input that cannot use it.
#
# Requires: Python 3 standard library, complex_variants.py alongside this file,
# and samtools on PATH when using --bam.
# ----------------------------------------------------------------------------

import argparse
import bisect
import gzip
import os
import re
import subprocess
import sys

import complex_variants as cv

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')
SIMPLE_CIGAR_RE = re.compile(r'^[0-9MX=]+$')

F_UNMAPPED = 0x4
F_SECONDARY = 0x100
F_SUPPLEMENTARY = 0x800

MASK = (1 << 64) - 1
_PRIMARY = set(["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY"])


# ---------------------------------------------------------------------------
# membership set
# ---------------------------------------------------------------------------

def _fnv(chrom, pos, ref, alt):
    """FNV-1a over the joined key. Deterministic across runs, unlike Python's
    salted hash() on str."""
    b = ("%s:%s:%s:%s" % (chrom, pos, ref, alt)).encode()
    x = 0xcbf29ce484222325
    for c in b:
        x ^= c
        x = (x * 0x100000001b3) & MASK
    return x


class RefFasta:
    """Random access to a .fai-indexed FASTA, holding one contig at a time.

    createDiff runs per contig and the BAM is coordinate-sorted, so one resident
    contig suffices; the 'other'/'unplaced' buckets hold many small contigs,
    which are cheap to swap. Reads the .fai directly, so no pysam dependency.
    """

    def __init__(self, path):
        self.path = path
        self.idx = {}
        with open(path + ".fai") as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                # name -> (length, offset, bases_per_line, bytes_per_line)
                self.idx[c[0]] = (int(c[1]), int(c[2]), int(c[3]), int(c[4]))
        self.fh = open(path, "rb")
        self.cur = None
        self.seq = None

    def _load(self, chrom):
        meta = self.idx.get(chrom)
        self.cur = chrom
        if meta is None:
            self.seq = None
            return
        length, offset, lb, lw = meta
        nlines = (length + lb - 1) // lb
        self.fh.seek(offset)
        raw = self.fh.read(nlines * lw)
        self.seq = raw.replace(b"\n", b"").replace(b"\r", b"")[:length].upper()

    def get(self, chrom, start0, length):
        """0-based half-open slice as str, or None if out of range/unknown."""
        if chrom != self.cur:
            self._load(chrom)
        if self.seq is None or start0 < 0 or length <= 0:
            return None
        if start0 + length > len(self.seq):
            return None
        return self.seq[start0:start0 + length].decode("ascii")


MAX_SHIFT = 500          # repeats longer than this are not worth chasing


def left_align(chrom, pos, ref, alt, ref_get):
    """Left-align one indel to match `bcftools norm -f` (the vt algorithm).

    pos is 1-based; returns (pos, ref, alt). Right-trim shared trailing bases,
    then roll left only while both alleles still end in the same base, prepending
    the reference base at pos-1. Substitutions/MNPs and anything without
    reference context are returned unchanged, which keeps the caller fail-closed.

    Validated against bcftools norm on real chr20 indels (see
    validate_norm_vs_bcftools.sh): must agree exactly.
    """
    if len(ref) == len(alt):
        return pos, ref, alt
    # right-trim shared trailing bases, keeping at least one base per allele
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    # roll left while both alleles end in the same base
    shifts = 0
    while ref[-1] == alt[-1] and pos > 1 and shifts < MAX_SHIFT:
        b = ref_get(chrom, pos - 2, 1)        # base just left of pos, 0-based
        if not b or b not in "ACGT":
            break
        ref = b + ref[:-1]
        alt = b + alt[:-1]
        pos -= 1
        shifts += 1
    return pos, ref, alt


def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def load_allowed(paths):
    """5-column BEDs (CHROM START0 END REF ALT) -> set of 64-bit hashes keyed on
    (chrom, END, ref, alt). END is the 1-based position, matching the payload's
    coordinates. Hashed rather than stored as strings: the union of the allowed
    list and the error set runs to ~2e8 alleles, which is ~22 GB as strings and
    ~7 GB as hashes."""
    allowed = set()
    for path in paths:
        if not path:
            continue
        with open_maybe_gz(path) as fh:
            for line in fh:
                if not line or line[0] == "#":
                    continue
                c = line.rstrip("\n").split("\t")
                if len(c) < 5:
                    continue
                allowed.add(_fnv(c[0], c[2], c[3], c[4]))
    return allowed


def load_blacklist(path):
    """{chrom: (starts[], ends[])}, sorted, for membership via bisect.
    BED is 0-based half-open."""
    if not path:
        return None
    iv = {}
    with open_maybe_gz(path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 3:
                continue
            iv.setdefault(c[0], []).append((int(c[1]), int(c[2])))
    bl = {}
    for chrom, lst in iv.items():
        lst.sort()
        bl[chrom] = ([a for a, _ in lst], [b for _, b in lst])
    return bl


def in_blacklist(bl, chrom, pos1):
    """pos1 is 1-based; BED is 0-based half-open."""
    if bl is None or chrom not in bl:
        return False
    starts, ends = bl[chrom]
    pos0 = pos1 - 1
    i = bisect.bisect_right(starts, pos0) - 1
    return i >= 0 and pos0 < ends[i]


def is_primary(chrom):
    return chrom in _PRIMARY


# ---------------------------------------------------------------------------
# substitution extraction
# ---------------------------------------------------------------------------

def parse_md(md):
    out = []
    i = 0
    n = len(md)
    while i < n:
        if md[i].isdigit():
            j = i
            while j < n and md[j].isdigit():
                j += 1
            out.append(("m", int(md[i:j])))
            i = j
        elif md[i] == "^":
            j = i + 1
            while j < n and md[j].isalpha():
                j += 1
            out.append(("d", md[i + 1:j]))
            i = j
        else:
            out.append(("x", md[i]))
            i += 1
    return out


def extract_subs(pos, seq, qual, md, min_bq):
    """Substitutions as (genomic_pos, ref_base, alt_base) for an M-only read.

    A mismatch whose read base is below min_bq is treated as a sequencing error
    and NOT recorded, so the read reconstructs as reference there. min_bq=0
    disables this; the v3 pipeline runs at 0 and identifies errors downstream
    by allele fraction, which base quality cannot do."""
    have_qual = qual != "*" and len(qual) == len(seq)
    subs = []
    ref_pos = pos
    read_idx = 0
    for kind, val in parse_md(md):
        if kind == "m":
            ref_pos += val
            read_idx += val
        elif kind == "d":
            ref_pos += len(val)
        else:
            if read_idx < len(seq):
                if not (min_bq > 0 and have_qual
                        and (ord(qual[read_idx]) - 33) < min_bq):
                    subs.append((ref_pos, val, seq[read_idx]))
            ref_pos += 1
            read_idx += 1
    return subs


def get_md(fields):
    for t in fields[11:]:
        if t.startswith("MD:Z:"):
            return t[5:]
    return None


# ---------------------------------------------------------------------------
# classify + route
# ---------------------------------------------------------------------------

def classify(fields, min_bq):
    """-> (cls, payload_cols, variants) where variants is a list of
    (pos, ref, alt) used for routing. cls None means reference-only."""
    flag = int(fields[1])
    cigar = fields[5]
    seq = fields[9]
    qual = fields[10]

    # Only AS is carried. This matches the previous format and is not an
    # oversight: pbam2bam emits recomputed MD/NM, the stored AS, and RG from
    # the header, and discards everything else. Storing the full tag set costs
    # ~100 bytes per record that nothing ever reads -- measured on chr2, that
    # was 210 B/record against 107 for the old format.
    #
    # The consequence is a real one and belongs in Methods: tags other than
    # MD/NM/AS/RG do NOT survive a round trip. "Recovers the original file
    # exactly" is true of the core fields and QUAL, not of the tag set.
    restore = []
    for t in fields[11:]:
        if t.startswith("AS:"):
            restore = [t]
            break
    mapped = not (flag & F_UNMAPPED)
    md = get_md(fields)

    # S: mapped, primary, M-only CIGAR, MD present. Substitutions are stored;
    # SEQ is not, since it rebuilds from the reference plus the substitutions.
    if (mapped and md is not None and seq != "*"
            and not (flag & (F_SECONDARY | F_SUPPLEMENTARY))
            and SIMPLE_CIGAR_RE.match(cigar)):
        subs = extract_subs(int(fields[3]), seq, qual, md, min_bq)
        if not subs:
            return None, None, None            # reference-only: no record
        sfield = ";".join("%d:%s>%s" % (p, r, a) for p, r, a in subs)
        return "S", ["S", sfield] + restore, subs

    # C: mapped with MD and real SEQ. Full edit set; SEQ and QUAL both omitted.
    # CIGAR IS stored: sanitize alters it, and it is the only such field.
    if mapped and md is not None and seq != "*":
        ex = cv.extract(int(fields[3]), cigar, seq, md, qual, min_bq)
        variants = ex["subs"] + ex["indels"]
        vfield = ";".join("%d:%s>%s" % (p, r, a) for p, r, a in variants) \
            if variants else "."
        splices = ["%d:I:%s" % (o, b) for o, b in ex["inserts"]] + \
                  ["%d:S:%s" % (o, b) for o, b in ex["clips"]]
        sfield = ";".join(splices) if splices else "."
        fmp = "1" if ex["force_mp"] else "0"
        return "C", ["C", cigar, vfield, sfield, fmp] + restore, variants

    # U: unmapped or no MD. Cannot rebuild from the reference, so store SEQ.
    return "U", ["U", seq] + restore, None


def route(cls, fields, variants, forcemp, allowed, bl, args):
    """-> True for MP, False for LP. Fail-closed: MP unless every difference is
    shown to be non-characterizing.

    Order matters and matches the previous implementation: the whole-read
    escapes are tested BEFORE forcemp, so a low-MAPQ or non-primary complex read
    goes to LP even when forcemp is set. Testing forcemp first would send it to
    MP instead, which is a different partition."""
    if cls == "U":
        return True                                   # unmapped: always MP

    chrom = fields[2]
    mapq = int(fields[4]) if fields[4].isdigit() else 0

    # Whole-read escapes. A read on a decoy/alt contig is not primary-assembly
    # genotype, and a read below the MAPQ floor is not confidently placed, so
    # neither carries characterizing information.
    if args.primary_chroms_only and not is_primary(chrom):
        return False
    if args.min_mapq > 0 and mapq < args.min_mapq:
        return False

    if forcemp:
        return True                    # soft clip / unanchorable indel

    rg = getattr(args, "_ref", None)
    for pos, ref, alt in (variants or []):
        if alt == "N":
            continue                                  # no-call: no genotype
        p, r, a = pos, ref, alt
        if rg is not None and len(ref) != len(alt):
            # Canonicalize before the membership test so the read's indel is
            # compared in the spelling the catalogue uses. Deliberately a
            # REPLACEMENT, not a fallback: testing the as-placed form first
            # would still let a rare indel whose placement collides with a
            # common one clear into LP, which is the leak this closes.
            p, r, a = left_align(chrom, int(pos), ref, alt, rg.get)
            if args.norm_audit is not None and (p, r, a) != (int(pos), ref, alt):
                args.norm_audit.write("%s\t%s\t%s>%s\t%s\t%s>%s\n"
                                      % (chrom, pos, ref, alt, p, r, a))
        if _fnv(chrom, p, r, a) in allowed:
            continue                                  # common, or a known error
        if in_blacklist(bl, chrom, int(pos)):
            continue                                  # unmappable region
        return True
    return False


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        prog="createDiff.py",
        description="Read a BAM and write the LP and MP layers, with records "
                    "keyed by pBAM line number.")
    src = ap.add_argument_group("input")
    src.add_argument("--bam", default=None,
                     help="input BAM (needs samtools). Default: SAM on stdin.")
    src.add_argument("--region", default=None,
                     help="restrict to a region, e.g. chr7. Line numbers then "
                          "run within the region; sanitize the SAME region so "
                          "the pBAM pairs with it.")
    src.add_argument("--min-bq", type=int, default=0,
                     help="Phred floor for calling a mismatch a substitution "
                          "(default 0 = off; errors are handled downstream by "
                          "allele fraction)")

    rt = ap.add_argument_group("routing")
    rt.add_argument("--allowed", default=None,
                    help="allowed-list BED(.gz): CHROM START0 END REF ALT. "
                         "Optional: without it nothing is non-characterizing, so "
                         "every variant-bearing read goes to MP. That is the "
                         "fail-closed default, and it skips a ~13 GB load for "
                         "input that cannot use it (e.g. unmapped-only streams, "
                         "where class U routes to MP without consulting the set).")
    rt.add_argument("--errset", default=None,
                    help="error-site BED(.gz); unioned with --allowed")
    rt.add_argument("--blacklist", default=None,
                    help="problem-region BED(.gz); substitutions inside never force MP")
    rt.add_argument("--min-mapq", type=int, default=0,
                    help="reads below this MAPQ do not force MP (0 = off)")
    rt.add_argument("--primary-chroms-only", action="store_true",
                    help="reads off chr1-22,X,Y do not force MP")
    rt.add_argument("--reference", default=None,
                    help="indexed FASTA (.fai beside it). Enables indel "
                         "left-alignment before the allowed-list test, so a "
                         "read's indel is compared in the same spelling the "
                         "catalogue uses. Without it, indels are tested exactly "
                         "as the aligner placed them.")
    rt.add_argument("--norm-audit", default=None, metavar="FILE",
                    help="log every indel whose spelling changed under "
                         "left-alignment: chrom, as-placed, normalized")

    out = ap.add_argument_group("output")
    out.add_argument("--out-prefix", required=True,
                     help="writes <prefix>.LP.diff and <prefix>.MP.diff")
    out.add_argument("--gzip", action="store_true")
    args = ap.parse_args()

    if args.reference:
        if not os.path.exists(args.reference + ".fai"):
            sys.stderr.write("ERROR: no .fai beside %s (run samtools faidx)\n"
                             % args.reference)
            return 1
        args._ref = RefFasta(args.reference)
        sys.stderr.write("indel normalization: ON (%s)\n" % args.reference)
    else:
        args._ref = None
        sys.stderr.write("indel normalization: off (no --reference)\n")
    args.norm_audit = open(args.norm_audit, "w") if args.norm_audit else None

    if args.allowed or args.errset:
        sys.stderr.write("loading membership set...\n")
        allowed = load_allowed([args.allowed, args.errset])
        sys.stderr.write("  alleles: %d\n" % len(allowed))
    else:
        allowed = set()
        sys.stderr.write("no --allowed/--errset: every variant-bearing read -> MP\n")
    bl = load_blacklist(args.blacklist)
    if bl is not None:
        sys.stderr.write("  blacklist chroms: %d\n" % len(bl))

    if args.bam:
        cmd = ["samtools", "view", args.bam]
        if args.region:
            cmd.append(args.region)
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        stream = proc.stdout
    else:
        if args.region:
            sys.stderr.write("--region requires --bam\n")
            return 1
        proc = None
        stream = sys.stdin

    suffix = ".gz" if args.gzip else ""
    opener = (lambda p: gzip.open(p, "wt")) if args.gzip else (lambda p: open(p, "wt"))
    lp_path = args.out_prefix + ".LP.diff" + suffix
    mp_path = args.out_prefix + ".MP.diff" + suffix

    lineno = 0
    last = {"LP": 0, "MP": 0}
    n = {"LP": 0, "MP": 0, "ref": 0}
    ncls = {"S": 0, "C": 0, "U": 0}

    with opener(lp_path) as lp, opener(mp_path) as mp:
        handle = {"LP": lp, "MP": mp}
        for line in stream:
            if line[0] == "@":
                continue
            lineno += 1                       # 1-based, matches the pBAM line
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue

            cls, cols, variants = classify(f, args.min_bq)
            if cls is None:
                n["ref"] += 1                 # reference-only: no record
                continue
            ncls[cls] += 1

            # C payload is ["C", cigar, variants, splices, forcemp], so cols[4]
            forcemp = (cls == "C" and cols[4] == "1")
            layer = "MP" if route(cls, f, variants, forcemp, allowed, bl, args) \
                    else "LP"

            # delta is per-file: the gap since the last record in THIS layer
            handle[layer].write("%d\t%s\n" % (lineno - last[layer],
                                              "\t".join(cols)))
            last[layer] = lineno
            n[layer] += 1

    if proc is not None:
        proc.stdout.close()
        proc.wait()

    if args.norm_audit is not None:
        args.norm_audit.close()

    sys.stderr.write(
        "done.\n"
        "  reads seen:    %d\n"
        "  reference-only: %d  (no record; reconstructed from the pBAM)\n"
        "  LP records:    %d\n"
        "  MP records:    %d\n"
        "  classes: S=%d C=%d U=%d\n"
        % (lineno, n["ref"], n["LP"], n["MP"], ncls["S"], ncls["C"], ncls["U"]))
    if lineno != n["ref"] + n["LP"] + n["MP"]:
        sys.stderr.write("  WARNING: reads seen != ref + LP + MP. Every read must\n"
                         "  land in exactly one bucket.\n")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
