#!/usr/bin/env python3
# ============================================================================
# partition_diff.py
#
# Part of the uni-id partitioning toolkit. Splits a keyed diff (produced by
# createDiff.py) into two privacy-stratified layers:
#
#   LP (less-private) : reads whose every variant is common, a sequencing
#                       error, a no-call, or otherwise non-characterizing.
#   MP (most-private) : reads carrying at least one rare / novel / identifying
#                       variant (the fail-closed default).
#
# A read is sent to LP only if EVERY substitution it carries is non-forcing.
# A substitution is non-forcing when any of the following hold:
#   * it is on the allowed list (common variant, or a supplied error site),
#   * its ALT base is N (a no-call: not a real, characterizing variant),
#   * it falls inside a supplied blacklist region (mismapping-prone),
#   * (optionally) the whole read is on a non-primary contig, or has MAPQ
#     below a threshold -- both signals of unreliable placement.
# Anything else forces the read to MP.
#
# This file provides two subcommands:
#   build-errset  Pile up one chromosome of a BAM and emit the sequencing-error
#                 sites (low-VAF, adequate-depth substitutions) as a BED. Run
#                 once per chromosome (parallelize across chromosomes), then
#                 concatenate the outputs into a single error-site BED.
#   partition     Split a diff into LP/MP against an allowed list (optionally
#                 unioned with an error-site BED, a blacklist, and MAPQ /
#                 primary-contig rules).
#
# The membership set is stored as 64-bit FNV-1a hashes rather than full
# "chrom:pos:ref:alt" strings, which cuts peak memory ~3x (a ~200M-entry set
# fits in ~7GB instead of ~22GB). A 64-bit collision over ~200M entries has
# probability ~2e-3 and can only ever route ONE read to LP that should have
# gone to MP -- harmless for the error set and astronomically rare otherwise.
#
# Requirements:
#   build-errset : bcftools (>=1.12) and samtools on PATH.
#   partition    : Python 3 standard library only.
# ============================================================================

import sys
import gzip
import bisect
import argparse
import subprocess


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

MASK = (1 << 64) - 1


def _fnv(chrom, pos, ref, alt):
    """FNV-1a over the joined key bytes; deterministic across runs/processes
    (unlike Python's salted hash() on str)."""
    b = ("%s:%s:%s:%s" % (chrom, pos, ref, alt)).encode()
    x = 0xcbf29ce484222325
    for c in b:
        x ^= c
        x = (x * 0x100000001b3) & MASK
    return x


def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


_PRIMARY = set(["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY"])


def is_primary(chrom):
    return chrom in _PRIMARY


# ---------------------------------------------------------------------------
# build-errset subcommand
# ---------------------------------------------------------------------------
# Pile up one chromosome and emit every substitution ALLELE that looks like a
# sequencing error: variant allele fraction (VAF) below --vaf-max at a site
# with depth at least --dp-min. Output is a 5-column BED
#   CHROM  START0  END  REF  ALT
# matching the allowed-list format, so the two can simply be concatenated into
# the membership set that `partition` consumes.
#
# The pileup is intentionally permissive so that its read set matches the diff:
#   -A            count anomalous (orphan) read pairs
#   --ff UNMAP    only exclude unmapped reads (keep duplicates, secondary, etc.)
#   -Q 0          no base-quality filtering (VAF, not quality, defines an error)
#   -d 10000      raise the depth cap (the DEFAULT caps low; note that -d 0 does
#                 NOT mean "no cap" -- it subsamples to ~zero, so never use it)
#   --no-BAQ      disable base-alignment-quality recomputation
#
# Only single-base REF and single-base ALT records are emitted (SNVs); the
# spanning-deletion marker "<*>" and multi-base alleles are skipped.

def run_build_errset(args):
    region = args.chrom
    mpileup = [
        "bcftools", "mpileup",
        "-f", args.reference,
        "-r", region,
        "-A", "--ff", "UNMAP",
        "-Q", "0",
        "-d", str(args.max_depth),
        "--no-BAQ",
        "-a", "FORMAT/AD,FORMAT/DP",
        args.bam,
    ]
    query = [
        "bcftools", "query",
        "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%AD]\t[%DP]\n",
    ]

    sys.stderr.write("[build-errset] %s : VAF<%g DP>=%d\n"
                     % (region, args.vaf_max, args.dp_min))

    p1 = subprocess.Popen(mpileup, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(query, stdin=p1.stdout, stdout=subprocess.PIPE, text=True)
    p1.stdout.close()

    out = open_maybe_gz(args.out, "wt") if args.out.endswith(".gz") else open(args.out, "w")
    n = 0
    try:
        for line in p2.stdout:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            chrom, pos, ref, alt_field, ad_field, dp_field = cols[:6]
            try:
                dp = int(dp_field)
            except ValueError:
                continue
            if dp < args.dp_min:
                continue
            alts = alt_field.split(",")
            ads = ad_field.split(",")
            # AD is [ref, alt1, alt2, ...]; align alt i with ads[i+1].
            for i, a in enumerate(alts):
                if a in ("<*>", ".", "") or len(a) != 1 or len(ref) != 1:
                    continue
                try:
                    altc = int(ads[i + 1])
                except (ValueError, IndexError):
                    continue
                if altc <= 0:
                    continue
                if altc / dp < args.vaf_max:
                    out.write("%s\t%d\t%s\t%s\t%s\n" % (chrom, int(pos) - 1, pos, ref, a))
                    n += 1
    finally:
        out.close()
        p2.stdout.close()
        p2.wait()
        p1.wait()

    sys.stderr.write("[build-errset] %s : %d error alleles -> %s\n" % (region, n, args.out))


# ---------------------------------------------------------------------------
# partition subcommand
# ---------------------------------------------------------------------------

def load_blacklist(path):
    """Return {chrom: (starts[], ends[])} sorted, for interval membership via
    bisect. BED is 0-based half-open [start, end)."""
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
            chrom, s, e = c[0], int(c[1]), int(c[2])
            iv.setdefault(chrom, []).append((s, e))
    bl = {}
    for chrom, lst in iv.items():
        lst.sort()
        bl[chrom] = ([a for a, _ in lst], [b for _, b in lst])
    return bl


def in_blacklist(bl, chrom, pos1):
    """pos1 is 1-based; BED is 0-based half-open. Convert to pos0 = pos1 - 1."""
    if bl is None or chrom not in bl:
        return False
    starts, ends = bl[chrom]
    pos0 = pos1 - 1
    i = bisect.bisect_right(starts, pos0) - 1
    return i >= 0 and pos0 < ends[i]


def load_allowed(paths):
    """Load one or more 5-column BEDs (CHROM START0 END REF ALT) into a set of
    64-bit hashes keyed on (chrom, END, ref, alt). END is the 1-based position,
    matching the diff payload's coordinate. Multiple paths are unioned (e.g. an
    allowed list plus an error-site BED)."""
    allowed = set()
    for path in paths:
        with open_maybe_gz(path) as fh:
            for line in fh:
                if not line or line[0] == "#":
                    continue
                c = line.rstrip("\n").split("\t")
                if len(c) < 5:
                    continue
                chrom, _s, end, ref, alt = c[0], c[1], c[2], c[3], c[4]
                allowed.add(_fnv(chrom, end, ref, alt))
    return allowed


def parse_subs(payload):
    if payload == ".":
        return []
    out = []
    for item in payload.split(";"):
        pos, ra = item.split(":")
        ref, alt = ra.split(">")
        out.append((pos, ref, alt))
    return out


def run_partition(args):
    allowed_paths = [args.allowed]
    if args.errset:
        allowed_paths.append(args.errset)

    sys.stderr.write("loading membership set (hashed) from %d file(s)...\n" % len(allowed_paths))
    allowed = load_allowed(allowed_paths)
    sys.stderr.write("membership alleles: %d\n" % len(allowed))
    bl = load_blacklist(args.blacklist)
    if bl is not None:
        sys.stderr.write("blacklist chroms: %d\n" % len(bl))

    suffix = ".gz" if args.gzip else ""
    lp_path = args.out_prefix + ".LP.diff" + suffix
    mp_path = args.out_prefix + ".MP.diff" + suffix
    opener = (lambda p: gzip.open(p, "wt")) if args.gzip else (lambda p: open(p, "wt"))

    n_lp = n_mp = n_c = n_ref = n_snv_lp = n_snv_mp = n_c_lp = 0

    with open_maybe_gz(args.diff) as din, opener(lp_path) as lp, opener(mp_path) as mp:
        for line in din:
            if line[0] == "#":
                lp.write(line); mp.write(line); continue
            f = line.rstrip("\n").split("\t")
            cls = f[9]; payload = f[10]

            # U (raw / unmapped): always MP.
            if cls == "U":
                mp.write(line); n_mp += 1; n_c += 1; continue

            # C (structured complex): route by its variant set, subject to the
            # same non-forcing rules as simple reads, plus a whole-read escape.
            if cls == "C":
                n_c += 1
                forcemp = (len(f) > 12 and f[12] == "1")
                chrom = f[2]
                read_mapq = int(f[4]) if f[4].isdigit() else 0
                if (args.primary_chroms_only and not is_primary(chrom)) or \
                   (args.min_mapq > 0 and read_mapq < args.min_mapq):
                    lp.write(line); n_lp += 1; n_c_lp += 1
                    continue
                variants = parse_subs(f[10])
                all_allowed = (not forcemp)
                if all_allowed:
                    for pos, ref, alt in variants:
                        if alt == "N":
                            continue
                        if _fnv(chrom, pos, ref, alt) in allowed:
                            continue
                        if in_blacklist(bl, chrom, int(pos)):
                            continue
                        all_allowed = False; break
                if all_allowed:
                    lp.write(line); n_lp += 1; n_c_lp += 1
                else:
                    mp.write(line); n_mp += 1
                continue

            # Reference-only simple read (no substitutions).
            if payload == ".":
                if args.keep_ref_only:
                    lp.write(line); n_lp += 1
                n_ref += 1
                continue

            # S / Sq (simple, with substitutions).
            chrom = f[2]
            read_mapq = int(f[4]) if f[4].isdigit() else 0
            # whole-read escapes: non-primary contig or sub-threshold MAPQ.
            if (args.primary_chroms_only and not is_primary(chrom)) or \
               (args.min_mapq > 0 and read_mapq < args.min_mapq):
                lp.write(line); n_lp += 1; n_snv_lp += 1
                continue
            subs = parse_subs(payload)
            all_allowed = True
            for pos, ref, alt in subs:
                if alt == "N":
                    continue  # no-call base: non-characterizing
                if _fnv(chrom, pos, ref, alt) in allowed:
                    continue  # common or confirmed-error: does not force MP
                if in_blacklist(bl, chrom, int(pos)):
                    continue  # problem region: mismapping artifact
                all_allowed = False; break
            if all_allowed:
                lp.write(line); n_lp += 1; n_snv_lp += 1
            else:
                mp.write(line); n_mp += 1; n_snv_mp += 1

    sys.stderr.write(
        "done.\n"
        "  ref-only reads: %d\n"
        "  LP records: %d (allowed-SNV %d, allowed-complex %d)\n"
        "  MP records: %d (complex seen %d)\n"
        % (n_ref, n_lp, n_snv_lp, n_c_lp, n_mp, n_c)
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser():
    ap = argparse.ArgumentParser(
        prog="partition_diff.py",
        description="Build a sequencing-error site set, or partition a keyed diff into LP/MP.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    # build-errset
    be = sub.add_parser("build-errset",
                        help="pile up one chromosome and emit low-VAF error sites as a BED")
    be.add_argument("--bam", required=True, help="input BAM/CRAM (indexed)")
    be.add_argument("--reference", required=True, help="reference FASTA (indexed)")
    be.add_argument("--chrom", required=True, help="chromosome / region to pile up, e.g. chr1")
    be.add_argument("--out", required=True, help="output BED(.gz): CHROM START0 END REF ALT")
    be.add_argument("--vaf-max", type=float, default=0.15,
                    help="a substitution with VAF below this is an error (default: 0.15)")
    be.add_argument("--dp-min", type=int, default=10,
                    help="minimum site depth to consider (default: 10)")
    be.add_argument("--max-depth", type=int, default=10000,
                    help="bcftools mpileup -d cap (default: 10000; do NOT use 0)")
    be.set_defaults(func=run_build_errset)

    # partition
    pa = sub.add_parser("partition", help="split a keyed diff into LP and MP layers")
    pa.add_argument("--diff", required=True, help="keyed diff from createDiff.py (.gz ok)")
    pa.add_argument("--allowed", required=True,
                    help="allowed-list BED(.gz): common variants, CHROM START0 END REF ALT")
    pa.add_argument("--errset", default=None,
                    help="optional error-site BED(.gz) (from build-errset); unioned with --allowed")
    pa.add_argument("--out-prefix", required=True, help="output prefix; writes <prefix>.LP.diff / .MP.diff")
    pa.add_argument("--gzip", action="store_true", help="gzip the output diffs")
    pa.add_argument("--blacklist", default=None,
                    help="BED(.gz) of problem regions; substitutions inside never force MP")
    pa.add_argument("--min-mapq", type=int, default=0,
                    help="reads with MAPQ below this do not force MP (0 = off)")
    pa.add_argument("--primary-chroms-only", action="store_true",
                    help="substitutions on non-primary contigs (not chr1..22,X,Y) never force MP")
    pa.add_argument("--keep-ref-only", action="store_true",
                    help="also emit reference-only simple reads into LP (default: drop them)")
    pa.set_defaults(func=run_partition)

    return ap


def main():
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
