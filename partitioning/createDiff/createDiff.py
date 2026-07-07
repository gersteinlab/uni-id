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
# Keyed, lossless replacement for ptools createDiff.
#
# Reads SAM records on stdin (samtools view <bam> | createDiff.py) and writes
# one tab-separated diff record per read on stdout. Unlike stock ptools, this captures
# the read's ALT bases explicitly, so substitutions round-trip on plain M-CIGAR reads
# (stock ptools loses them — it reconstructs the reference base at every mismatch).
#
# Two record classes:
#   S (simple)  : mapped, primary, CIGAR is only M/=/X, and an MD tag is present.
#                 Substitutions are extracted as genomic (pos:ref>alt) so the partition
#                 step can test them directly against allowed_loci.bed.gz with no offset
#                 math. SEQ is NOT stored — it is rebuilt at reconstruction from the
#                 reference plus these substitutions.
#   C (complex) : everything else (indels, soft/hard clips, N, unmapped, secondary/
#                 supplementary, missing MD). These are routed to MP unconditionally in
#                 v1. The full original SEQ and QUAL are stored so they round-trip
#                 without any indel-reconstruction logic (deferred to v2).
#
# Output columns (tab-separated; RESTORE is the final field and may itself contain tabs):
#   0  QNAME
#   1  FLAG
#   2  RNAME
#   3  POS
#   4  MAPQ
#   5  CIGAR
#   6  RNEXT
#   7  PNEXT
#   8  TLEN
#   9  CLASS     S or C
#   10 PAYLOAD   S: "pos:ref>alt;..." (or "." if the read matches reference)
#                C: the raw original SEQ
#   11 QUAL      C: raw original QUAL ; S: "." (QUAL is taken from the pBAM on rebuild)
#   12+ RESTORE  the original optional tag columns (12.. of the SAM line), verbatim
#
# Columns 0-8 are the reconstruction key. They are invariant under the pBAM transform
# for simple M reads (the LP-relevant case). Complex reads carry their full sequence so
# they never depend on key-matching against an altered pBAM record.

import sys
import os
import re
import argparse
import subprocess
import tempfile
import complex_variants as cv

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')

# SAM FLAG bits
F_UNMAPPED = 0x4
F_SECONDARY = 0x100
F_SUPPLEMENTARY = 0x800


def cigar_ops(cigar):
    """Return list of (length, op) tuples; [] if cigar is '*' or malformed."""
    if cigar == "*" or cigar == "":
        return []
    return [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]


def cigar_is_simple(ops):
    """True iff every op is an aligned-match op (M/=/X) — no indels, clips, skips."""
    return len(ops) > 0 and all(op in ("M", "=", "X") for _, op in ops)


def find_md(tags):
    """Return the MD value (string after 'MD:Z:') or None."""
    for t in tags:
        if t.startswith("MD:Z:"):
            return t[5:]
    return None


# MD tokens: a run-length integer, a deletion (^ + bases), or a single mismatch ref base
MD_RE = re.compile(r'(\d+)|(\^[A-Za-z]+)|([A-Za-z])')


def extract_subs(pos, seq, qual, md, min_bq):
    """
    Walk an MD string over a pure M/=/X read and return a list of
    (genomic_pos, ref_base, alt_base) substitutions, or None if the walk is
    inconsistent (caller then treats the read as complex).

    Valid only when CIGAR has no I/D/S/H/N: read index and reference index advance
    together, so an MD match-run of n advances both by n, and an MD mismatch base
    consumes one base of each. A '^' deletion token means the CIGAR was not actually
    simple (inconsistent input) -> bail to complex.

    Base-quality filter: a mismatch whose read-base Phred quality is below min_bq is
    treated as NOT a confirmed variant and is skipped (the position is left as
    reference). This drops sequencing errors, which present as low-quality singleton
    mismatches, so they do not route an otherwise-common read to MP. min_bq <= 0
    disables the filter. QUAL '*' (absent) disables it for that read.
    """
    have_qual = qual != "*" and len(qual) == len(seq)
    read_idx = 0
    ref_pos = pos  # 1-based genomic position of the next aligned base
    subs = []
    for num, deletion, base in MD_RE.findall(md):
        if num:
            n = int(num)
            read_idx += n
            ref_pos += n
        elif deletion:
            return None  # deletion in MD but CIGAR was simple: inconsistent
        elif base:
            if read_idx >= len(seq):
                return None  # ran off the read: inconsistent
            # skip low-quality mismatches: treat as reference, not a variant
            if min_bq > 0 and have_qual and (ord(qual[read_idx]) - 33) < min_bq:
                read_idx += 1
                ref_pos += 1
                continue
            alt = seq[read_idx]
            subs.append((ref_pos, base.upper(), alt.upper()))
            read_idx += 1
            ref_pos += 1
    if read_idx != len(seq):
        return None  # MD did not span the whole read: inconsistent
    return subs


def is_simple_fields(f, min_bq):
    """True iff this SAM record is a primary, mapped, pure-M read with MD and real SEQ
    AND its substitutions extract cleanly. Mirrors make_record's 'simple' decision so the
    collision pre-pass and the emit pass agree on which reads are simple."""
    flag = int(f[1])
    rname = f[2]
    seq = f[9]
    if flag & (F_SECONDARY | F_SUPPLEMENTARY):
        return False
    if (flag & F_UNMAPPED) or rname == "*":
        return False
    if seq == "*":
        return False
    if not cigar_is_simple(cigar_ops(f[5])):
        return False
    md = find_md(f[11:])
    if md is None:
        return False
    return extract_subs(int(f[3]), seq, f[10], md, min_bq) is not None


def position_key(f):
    """The 8-field key minus QNAME: FLAG,RNAME,POS,MAPQ,CIGAR,RNEXT,PNEXT,TLEN.
    Invariant under the pBAM transform for simple reads, so it identifies a simple read's
    pBAM image without storing QNAME (recovered from the matched pBAM read)."""
    return "\t".join(f[1:9])


def compute_collisions(bam, min_bq, tmpdir):
    """
    One streaming pass over the BAM: emit the position-key of every simple read to a temp
    file, then sort|uniq -d (disk-backed, memory bounded) to find position-keys shared by
    two or more simple reads. Those reads must keep QNAME (cannot be disambiguated by
    position alone). Returns a set of colliding position-key strings.
    """
    keyfile = tempfile.NamedTemporaryFile(mode="w", dir=tmpdir, prefix="poskeys.",
                                          suffix=".tsv", delete=False)
    try:
        view = subprocess.Popen(["samtools", "view", bam], stdout=subprocess.PIPE, text=True)
        n_simple = 0
        for line in view.stdout:
            f = line.rstrip("\n").split("\t")
            if is_simple_fields(f, min_bq):
                keyfile.write(position_key(f) + "\n")
                n_simple += 1
        view.stdout.close()
        view.wait()
        keyfile.close()

        collisions = set()
        # sort | uniq -d  via subprocess; read the duplicate keys back
        sort = subprocess.Popen(["sort", keyfile.name], stdout=subprocess.PIPE, text=True,
                                env=dict(os.environ, LC_ALL="C"))
        uniq = subprocess.Popen(["uniq", "-d"], stdin=sort.stdout, stdout=subprocess.PIPE,
                                text=True)
        sort.stdout.close()
        for line in uniq.stdout:
            collisions.add(line.rstrip("\n"))
        uniq.wait()
        sort.wait()
        sys.stderr.write("collision pre-pass: %d simple reads, %d colliding position-keys "
                         "(those reads keep QNAME)\n" % (n_simple, len(collisions)))
        return collisions
    finally:
        try:
            os.unlink(keyfile.name)
        except OSError:
            pass


def make_record(line, min_bq, collisions=None):
    line = line.rstrip("\n")
    f = line.split("\t")
    qname, flag_s, rname, pos_s, mapq, cigar, rnext, pnext, tlen = f[:9]
    seq = f[9]
    qual = f[10]
    tags = f[11:]
    flag = int(flag_s)

    key = [qname, flag_s, rname, pos_s, mapq, cigar, rnext, pnext, tlen]
    # v2: store only AS (genuinely unrecomputable). MD/NM are recomputed at reconstruction;
    # RG is re-stamped from the header; all other tags (PG, MQ, XS, MC, SA, XA, ...) are
    # dropped. This is the bulk of the v2 size win and is lossy on those tags by design.
    as_tag = next((t for t in tags if t.startswith("AS:")), None)
    restore = [as_tag] if as_tag else []

    # Decide class. Anything not cleanly a primary, mapped, pure-M read with an MD tag
    # and a real sequence is complex -> store raw SEQ/QUAL, route to MP.
    is_primary = not (flag & (F_SECONDARY | F_SUPPLEMENTARY))
    is_mapped = not (flag & F_UNMAPPED) and rname != "*"
    ops = cigar_ops(cigar)
    md = find_md(tags)

    simple = (
        is_primary
        and is_mapped
        and seq != "*"
        and cigar_is_simple(ops)
        and md is not None
    )

    if simple:
        subs = extract_subs(int(pos_s), seq, qual, md, min_bq)
        if subs is not None:
            payload = ";".join("%d:%s>%s" % (p, r, a) for p, r, a in subs) if subs else "."
            # QNAME-drop: if collisions is provided and this read's position-key is unique,
            # store class 'S' with QNAME blanked ('.'); reconstruction recovers QNAME from
            # the uniquely-matching pBAM read. If the position-key collides, or we have no
            # collision info (stdin mode), keep QNAME with class 'Sq'.
            if collisions is not None and position_key(f) not in collisions:
                return "\t".join(["."] + key[1:] + ["S", payload, "."] + restore)
            return "\t".join(key + ["Sq", payload, "."] + restore)
        # extraction was inconsistent: fall through to complex

    # complex. Two sub-cases:
    #   C (structured): mapped, has MD, real SEQ. Extract the full variant set (subs +
    #     indels) for classification and reconstruction. SEQ is NOT stored (rebuilt from
    #     reference at reconstruction); QUAL is NOT stored (recovered from the pBAM image).
    #       col 10  variants : "gpos:ref>alt;..." (subs AND indels; SNVs are ref/alt len 1)
    #       col 11  splices  : "offset:I:bases;offset:S:bases;..." (inserted/clipped bases)
    #       col 12  forcemp  : "1" if the read must stay MP (soft clip / unanchorable indel)
    #       col 13+ AS
    #   U (raw): unmapped or missing MD -- cannot rebuild from reference, so store raw SEQ
    #     (QUAL still recovered from the pBAM). col 10 = SEQ, col 11 = ".".
    if is_mapped and md is not None and seq != "*":
        ex = cv.extract(int(pos_s), cigar, seq, md, qual, min_bq)
        variants = ex["subs"] + ex["indels"]
        vfield = ";".join("%d:%s>%s" % (p, r, a) for p, r, a in variants) if variants else "."
        splices = [("%d:I:%s" % (o, b)) for o, b in ex["inserts"]] + \
                  [("%d:S:%s" % (o, b)) for o, b in ex["clips"]]
        sfield = ";".join(splices) if splices else "."
        fmp = "1" if ex["force_mp"] else "0"
        return "\t".join(key + ["C", vfield, sfield, fmp] + restore)

    # U: raw fallback (unmapped / no MD)
    return "\t".join(key + ["U", seq, "."] + restore)


def main():
    ap = argparse.ArgumentParser(
        description="Keyed, lossless createDiff with base-quality mismatch filtering and "
                    "QNAME-drop. With --bam it does a two-pass run (collision pre-pass + "
                    "emit) so simple reads with a unique position-key can drop QNAME. "
                    "Without --bam it reads SAM on stdin and keeps QNAME on every record."
    )
    ap.add_argument("--bam", help="input BAM/CRAM path; enables QNAME-drop (two-pass).")
    ap.add_argument("--min-bq", type=int, default=20,
                    help="Phred base-quality floor for calling a mismatch a substitution. "
                         "Mismatches below this are treated as reference. 0 disables. Default: 20")
    ap.add_argument("--tmpdir", default=".",
                    help="directory for the collision-sort temp file (default: current dir).")
    args = ap.parse_args()

    out = sys.stdout
    out.write("#QNAME\tFLAG\tRNAME\tPOS\tMAPQ\tCIGAR\tRNEXT\tPNEXT\tTLEN\tCLASS\tPAYLOAD\tQUAL\tRESTORE...\n")

    if args.bam:
        collisions = compute_collisions(args.bam, args.min_bq, args.tmpdir)
        view = subprocess.Popen(["samtools", "view", args.bam], stdout=subprocess.PIPE, text=True)
        for line in view.stdout:
            if not line or line[0] == "@":
                continue
            out.write(make_record(line, args.min_bq, collisions) + "\n")
        view.stdout.close()
        view.wait()
    else:
        # stdin filter mode: no collision info, QNAME kept on every record
        for line in sys.stdin:
            if not line or line[0] == "@":
                continue
            out.write(make_record(line, args.min_bq, None) + "\n")


if __name__ == "__main__":
    main()
