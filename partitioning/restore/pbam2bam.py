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
# pbam2bam.py
# Reconstruct a BAM (as SAM on stdout) from a pBAM plus one or more keyed diff layers
# (LP and/or MP) produced by createDiff.py + partition_diff.py.
#
# Access tiers fall out of which layers you pass:
#   pBAM only            -> public skeleton (all sequence is reference)
#   pBAM + LP            -> common-variant reads restored; rare/novel/tainted reads and
#                           all complex reads remain reference-sanitized (present-but-
#                           sanitized: the read is there, its sequence is reference)
#   pBAM + LP + MP       -> full original restored
#
# Mechanism (does NOT line-zip like stock ptools):
#   * Simple (S) records are loaded into a dict keyed by columns 0-8. These key fields
#     are invariant under the pBAM transform for pure-M reads, so each pBAM read is
#     matched by an exact key lookup.
#       - hit  : restore SEQ (reference + stored ALT bases at genomic_pos - POS offsets),
#                restore original tags, restore POS/CIGAR from the key.
#       - miss : emit the pBAM read AS-IS (sanitized). A miss is never an error; it means
#                the read's restoration lives in a layer not loaded, or it is the
#                sanitized image of a complex read (see below).
#   * Complex (C) records store the full original read (SEQ/QUAL/tags) and are emitted
#     DIRECTLY from the diff after the pBAM stream — never matched through the pBAM,
#     because the transform alters their POS/CIGAR. To avoid double-emitting, the
#     sanitized pBAM image of a held complex read is suppressed during the stream using
#     an invariant (QNAME, FLAG) sub-key. If a complex read's layer is NOT loaded, no C
#     record exists, nothing is suppressed, and its sanitized pBAM image passes through
#     (present-but-sanitized) -- exactly the LP-only privacy behavior.
#
# The pBAM is read on stdin as SAM (samtools view -h ... | pbam2bam.py ...), header
# lines (@) are passed through unchanged. Output is SAM on stdout.
#
# Usage:
#   samtools view -h pbam.bam \
#     | pbam2bam.py --diff test.LP.diff [--diff test.MP.diff] \
#     | samtools view -h -b - > reconstructed.bam

import sys
import re
import argparse
import complex_variants as cv

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')
NKEY = 9  # columns 0..8 are the reconstruction key


def load_layers(paths):
    """
    Read the diff layer files. Returns:
      simple_pk : dict position_key_tuple(f[1:9]) -> (payload, restore)   [class S, QNAME dropped]
      simple_fk : dict full_key_tuple(f[0:9])      -> (payload, restore)   [class Sq, QNAME kept]
      complex   : list of full C/U field-lists (emitted directly)
      csuppress : set of (QNAME, FLAG) for held complex reads (suppress their pBAM images)
      cqual_needed : invariant keys (QNAME,FLAG,RNEXT,PNEXT,TLEN) of complex reads
    """
    simple_pk = {}
    simple_fk = {}
    complex_records = []
    csuppress = set()
    cqual_needed = set()
    for path in paths:
        with open(path, "rt") as fh:
            for line in fh:
                if not line or line[0] == "#":
                    continue
                f = line.rstrip("\n").split("\t")
                cls = f[9]
                if cls == "S":          # QNAME dropped: key by position (f[1:9])
                    simple_pk[tuple(f[1:NKEY])] = (f[10], f[12:])
                elif cls == "Sq":       # QNAME kept: key by full 9-tuple
                    simple_fk[tuple(f[:NKEY])] = (f[10], f[12:])
                else:                   # C or U
                    complex_records.append(f)
                    csuppress.add((f[0], f[1]))
                    cqual_needed.add((f[0], f[1], f[6], f[7], f[8]))
    return simple_pk, simple_fk, complex_records, csuppress, cqual_needed


def apply_subs(ref_seq, pos, payload):
    """
    Given the reference-derived SEQ (from the pBAM read) and a substitution payload
    'gpos:ref>alt;...', write each ALT base back at offset (gpos - pos). Returns the
    restored SEQ. payload '.' -> unchanged.
    """
    if payload == ".":
        return ref_seq
    seq = list(ref_seq)
    for item in payload.split(";"):
        loc, ra = item.split(":")
        _ref, alt = ra.split(">")
        off = int(loc) - pos
        if 0 <= off < len(seq):
            seq[off] = alt
    return "".join(seq)


def rebuild_simple(pbam_fields, qname, rest8, payload, restore, rg_tag):
    """
    Rebuild a simple read's SAM line. qname is supplied separately (recovered from the
    pBAM read for class-S records, or the stored QNAME for class-Sq). rest8 is the other
    eight key fields (FLAG,RNAME,POS,MAPQ,CIGAR,RNEXT,PNEXT,TLEN). SEQ is the pBAM
    reference sequence with ALT bases written back; QUAL comes from the pBAM read; MD/NM
    are recomputed; tags are the restored originals.
    """
    flag, rname, pos_s, mapq, cigar, rnext, pnext, tlen = rest8
    ref_seq = pbam_fields[9]
    qual = pbam_fields[10]
    seq = apply_subs(ref_seq, int(pos_s), payload)
    md_tag, nm_tag = build_md_nm(payload, int(pos_s), len(seq))
    out = [qname, flag, rname, pos_s, mapq, cigar, rnext, pnext, tlen, seq, qual]
    out.append(md_tag)
    out.append(nm_tag)
    as_tag = as_from_restore(restore)
    if as_tag:
        out.append(as_tag)
    if rg_tag:
        out.append(rg_tag)
    return "\t".join(out)


def parse_splices(sfield):
    """'offset:I:bases;offset:S:bases' -> (inserts list, clips list) of (offset, bases)."""
    inserts, clips = [], []
    if sfield == "." or sfield == "":
        return inserts, clips
    for item in sfield.split(";"):
        off_s, typ, bases = item.split(":")
        off = int(off_s)
        if typ == "I":
            inserts.append((off, bases))
        elif typ == "S":
            clips.append((off, bases))
    return inserts, clips


def emit_complex(fields, qual, rg_tag, ref_lookup):
    """
    Emit a complex read with recovered QUAL and re-stamped RG.

    Class 'C' (structured): rebuild SEQ from reference + stored variants/splices.
      layout: 0..8 key, 9 'C', 10 variants(subs+indels), 11 splices, 12 forcemp, 13.. tags
    Class 'U' (raw): unmapped / no-MD; SEQ stored verbatim.
      layout: 0..8 key, 9 'U', 10 SEQ, 11 '.', 12.. tags

    QUAL is recovered from the pBAM image (col passed in). MD/NM are omitted for complex
    reads (acceptable lossy-tag policy). AS is restored if present.
    """
    cls = fields[9]
    key = fields[:NKEY]

    if cls == "U":
        seq = fields[10]
        tags = fields[12:]
    else:  # 'C'
        chrom = fields[2]
        pos = int(fields[3])
        cigar = fields[5]
        variants_field = fields[10]
        inserts, clips = parse_splices(fields[11])
        # SNV subs for M-position overwrite: entries with single-base ref AND alt
        subs = []
        if variants_field != ".":
            for it in variants_field.split(";"):
                loc, ra = it.split(":")
                ref, alt = ra.split(">")
                if len(ref) == 1 and len(alt) == 1:
                    subs.append((int(loc), ref, alt))
        ex = {"subs": subs, "inserts": inserts, "clips": clips}
        seq = cv.rebuild_seq(pos, cigar, ref_lookup, chrom, ex)
        tags = fields[13:]

    out = key + [seq, qual]
    as_tag = as_from_restore(tags)
    if as_tag:
        out.append(as_tag)
    if rg_tag:
        out.append(rg_tag)
    return "\t".join(out)


def build_md_nm(payload, pos, seqlen):
    """
    Recompute MD and NM for a simple (pure-M) read from its substitution payload.
    The subs ARE the complete set of mismatches vs reference, so this reproduces the
    original MD/NM exactly (lossless) without needing the reference FASTA.
    Returns (md_tag, nm_tag) as full SAM tag strings.
    """
    if payload == ".":
        return ("MD:Z:%d" % seqlen, "NM:i:0")
    items = []
    for it in payload.split(";"):
        loc, ra = it.split(":")
        ref, _alt = ra.split(">")
        items.append((int(loc) - pos, ref))
    items.sort()
    md = ""
    prev = 0
    for off, ref in items:
        md += str(off - prev) + ref
        prev = off + 1
    md += str(seqlen - prev)
    return ("MD:Z:" + md, "NM:i:%d" % len(items))


def parse_rg_id(header_line):
    """Extract the RG ID from an @RG header line, or None."""
    if not header_line.startswith("@RG"):
        return None
    for fld in header_line.rstrip("\n").split("\t"):
        if fld.startswith("ID:"):
            return fld[3:]
    return None


def as_from_restore(restore):
    """Return the stored AS tag (the only tag v2 keeps) or None."""
    for t in restore:
        if t.startswith("AS:"):
            return t
    return None


def main():
    ap = argparse.ArgumentParser(description="Reconstruct BAM from pBAM + keyed diff layers.")
    ap.add_argument("--diff", action="append", required=True,
                    help="a diff layer file (.LP.diff and/or .MP.diff). Repeatable.")
    ap.add_argument("--reference",
                    help="reference FASTA, required if any 'C' (structured complex) records "
                         "are present; their SEQ is rebuilt from it. 'U' (raw) and simple "
                         "reads do not need it.")
    args = ap.parse_args()

    simple_pk, simple_fk, complex_records, csuppress, cqual_needed = load_layers(args.diff)
    sys.stderr.write("loaded: %d simple(QNAME-dropped) + %d simple(QNAME-kept), %d complex "
                     "(suppress %d QNAME/FLAG)\n"
                     % (len(simple_pk), len(simple_fk), len(complex_records), len(csuppress)))

    # Load reference only if we actually have structured complex records to rebuild.
    n_struct = sum(1 for fl in complex_records if fl[9] == "C")
    ref_lookup = None
    if n_struct > 0:
        if not args.reference:
            sys.stderr.write("ERROR: %d structured complex (C) records require --reference\n"
                             % n_struct)
            sys.exit(2)
        sys.stderr.write("loading reference %s for %d complex reads...\n"
                         % (args.reference, n_struct))
        _ref = {}
        _name, _parts = None, []
        with open(args.reference) as fh:
            for line in fh:
                if line.startswith(">"):
                    if _name is not None:
                        _ref[_name] = "".join(_parts)
                    _name = line[1:].split()[0]
                    _parts = []
                else:
                    _parts.append(line.strip())
        if _name is not None:
            _ref[_name] = "".join(_parts)

        def ref_lookup(chrom, start0, length):
            return _ref[chrom][start0:start0 + length]

    out = sys.stdout
    n_hit = n_miss = n_suppressed = 0
    rg_tag = None  # first @RG ID seen in the pBAM header, re-stamped on every read
    complex_qual = {}  # invariant key -> QUAL, harvested from the pBAM stream

    for line in sys.stdin:
        if line[0] == "@":
            out.write(line)  # pass through SAM header
            if rg_tag is None:
                rg_id = parse_rg_id(line)
                if rg_id is not None:
                    rg_tag = "RG:Z:" + rg_id
            continue
        f = line.rstrip("\n").split("\t")

        # harvest QUAL for complex reads BEFORE any branch: the pBAM image of a complex
        # read is the one that gets suppressed, so grab its (verbatim) QUAL now.
        inv = (f[0], f[1], f[6], f[7], f[8])
        if inv in cqual_needed:
            complex_qual[inv] = f[10]

        # Suppress complex images FIRST, before the simple lookups. A complex read's pBAM
        # image has an altered position-key that could coincidentally collide with a real
        # simple read's position-key; checking (QNAME,FLAG) suppression first prevents it
        # from being mis-reconstructed as a simple read.
        if (f[0], f[1]) in csuppress:
            n_suppressed += 1
            continue

        # Sq (QNAME kept): match by full 9-key.
        rec = simple_fk.get(tuple(f[:NKEY]))
        if rec is not None:
            payload, restore = rec
            out.write(rebuild_simple(f, f[0], f[1:NKEY], payload, restore, rg_tag) + "\n")
            n_hit += 1
            continue

        # S (QNAME dropped): match by position-key (f[1:9]); recover QNAME from this pBAM read.
        rec = simple_pk.get(tuple(f[1:NKEY]))
        if rec is not None:
            payload, restore = rec
            out.write(rebuild_simple(f, f[0], f[1:NKEY], payload, restore, rg_tag) + "\n")
            n_hit += 1
            continue

        out.write(line)  # present-but-sanitized pass-through
        n_miss += 1

    # emit held complex reads directly, recovering QUAL by invariant key
    n_qual_missing = 0
    for fields in complex_records:
        inv = (fields[0], fields[1], fields[6], fields[7], fields[8])
        qual = complex_qual.get(inv)
        if qual is None:
            qual = "*"  # no pBAM image matched; emit unknown QUAL rather than fail
            n_qual_missing += 1
        out.write(emit_complex(fields, qual, rg_tag, ref_lookup) + "\n")

    sys.stderr.write(
        "done. pBAM reads: %d restored, %d passed-through-sanitized, %d suppressed; "
        "%d complex emitted from diff (%d with QUAL recovered, %d QUAL-missing).\n"
        % (n_hit, n_miss, n_suppressed, len(complex_records),
           len(complex_records) - n_qual_missing, n_qual_missing)
    )


if __name__ == "__main__":
    main()
