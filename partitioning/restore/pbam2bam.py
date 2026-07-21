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
import gzip
import argparse
import complex_variants as cv

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')
NKEY = 9  # columns 0..8 are the reconstruction key



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



# ---------------------------------------------------------------------------
# diff reader
# ---------------------------------------------------------------------------

def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


class DiffReader:
    """Streams a delta-encoded diff, exposing the absolute pBAM line of the next
    record. Records are read in order and never held: this is what makes the
    restore O(1) in memory instead of the 105-158 GB the old dict needed for
    175M LP records."""

    def __init__(self, path):
        self.fh = open_maybe_gz(path)
        self.path = path
        self.line = None      # absolute pBAM line of the pending record
        self.rec = None       # record columns after the delta
        self.n = 0
        self._abs = 0
        self.advance()

    def advance(self):
        for raw in self.fh:
            if not raw or raw[0] == "#":
                continue
            f = raw.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            self._abs += int(f[0])       # delta since this file's last record
            self.line = self._abs
            self.rec = f[1:]
            self.n += 1
            return
        self.line = None
        self.rec = None

    def close(self):
        self.fh.close()


# ---------------------------------------------------------------------------
# restore, one record class each
# ---------------------------------------------------------------------------
# Every one of these takes the pBAM read's fields for the parts sanitize
# preserved (QNAME, FLAG, RNAME, POS, MAPQ, RNEXT, PNEXT, TLEN, QUAL) and the
# diff record for the parts it did not.

def restore_simple(pf, rec, rg_tag):
    """rec = [S, subs, *tags]. The pBAM SEQ is the reference at this position,
    so restoring is writing the ALT bases back into it."""
    pos_s = pf[3]
    payload = rec[1]
    seq = apply_subs(pf[9], int(pos_s), payload)
    md_tag, nm_tag = build_md_nm(payload, int(pos_s), len(seq))
    out = pf[:9] + [seq, pf[10], md_tag, nm_tag]
    as_tag = as_from_restore(rec[2:])
    if as_tag:
        out.append(as_tag)
    if rg_tag:
        out.append(rg_tag)
    return "\t".join(out)


def restore_complex(pf, rec, ref_lookup, rg_tag):
    """rec = [C, cigar, variants, splices, forcemp, *tags].

    CIGAR comes from the record: it is the only field sanitize alters, which is
    measured (200/200 chr6 complex reads differ on CIGAR, 0/200 differ on
    anything else). Everything else comes from the pBAM line."""
    chrom = pf[2]
    pos = int(pf[3])
    cigar = rec[1]
    inserts, clips = parse_splices(rec[3])
    subs = []
    if rec[2] != ".":
        for it in rec[2].split(";"):
            loc, ra = it.split(":")
            ref, alt = ra.split(">")
            if len(ref) == 1 and len(alt) == 1:   # SNVs only; indels ride the CIGAR
                subs.append((int(loc), ref, alt))
    ex = {"subs": subs, "inserts": inserts, "clips": clips}
    seq = cv.rebuild_seq(pos, cigar, ref_lookup, chrom, ex)
    key = list(pf[:9])
    key[5] = cigar                                 # swap the sanitized CIGAR out
    out = key + [seq, pf[10]]
    as_tag = as_from_restore(rec[5:])
    if as_tag:
        out.append(as_tag)
    if rg_tag:
        out.append(rg_tag)
    return "\t".join(out)


def restore_unmapped(pf, rec, rg_tag):
    """rec = [U, seq, *tags]. Not rebuildable from the reference, so SEQ was
    stored verbatim. QUAL still comes from the pBAM."""
    out = pf[:9] + [rec[1], pf[10]]
    as_tag = as_from_restore(rec[2:])
    if as_tag:
        out.append(as_tag)
    if rg_tag:
        out.append(rg_tag)
    return "\t".join(out)


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        prog="pbam2bam.py",
        description="Restore a BAM from a pBAM plus one or more line-numbered "
                    "diff layers. Reads pBAM SAM on stdin, writes SAM on stdout.")
    ap.add_argument("--diff", action="append", required=True,
                    help="a layer (repeatable). pBAM+LP restores common variation; "
                         "pBAM+LP+MP restores the original.")
    ap.add_argument("--reference", default=None,
                    help="reference FASTA. Required only if a supplied layer holds "
                         "class-C records.")
    args = ap.parse_args()

    readers = [DiffReader(p) for p in args.diff]
    for r in readers:
        sys.stderr.write("layer %s: first record at pBAM line %s\n"
                         % (r.path, r.line))

    ref_lookup = None
    if args.reference:
        sys.stderr.write("loading reference %s...\n" % args.reference)
        _ref = {}
        _name, _parts = None, []
        opener = gzip.open if args.reference.endswith(".gz") else open
        with opener(args.reference, "rt") as fh:
            for l in fh:
                if l.startswith(">"):
                    if _name is not None:
                        _ref[_name] = "".join(_parts)
                    _name = l[1:].split()[0]
                    _parts = []
                else:
                    _parts.append(l.strip())
        if _name is not None:
            _ref[_name] = "".join(_parts)
        sys.stderr.write("  contigs: %d\n" % len(_ref))

        def ref_lookup(chrom, start0, length):
            return _ref[chrom][start0:start0 + length]

    n_restored = n_through = 0
    n_cls = {"S": 0, "C": 0, "U": 0}
    lineno = 0
    rg_tag = None      # first @RG ID in the pBAM header, re-stamped on every read

    for line in sys.stdin:
        if line[0] == "@":
            sys.stdout.write(line)
            if rg_tag is None:
                rg_id = parse_rg_id(line)
                if rg_id is not None:
                    rg_tag = "RG:Z:" + rg_id
            continue

        lineno += 1
        pf = line.rstrip("\n").split("\t")

        # Records are disjoint across layers by construction, so at most one
        # reader can be pending at this line. No key, so no collision is
        # representable: this is the whole point of the format.
        hit = None
        for r in readers:
            if r.line == lineno:
                hit = r
                break

        if hit is None:
            sys.stdout.write(line)          # reference-only: pass through
            n_through += 1
            continue

        rec = hit.rec
        hit.advance()
        cls = rec[0]
        n_cls[cls] = n_cls.get(cls, 0) + 1

        if cls == "S":
            sys.stdout.write(restore_simple(pf, rec, rg_tag) + "\n")
        elif cls == "U":
            sys.stdout.write(restore_unmapped(pf, rec, rg_tag) + "\n")
        else:
            if ref_lookup is None:
                sys.stderr.write("class-C record at line %d but no --reference\n"
                                 % lineno)
                return 1
            sys.stdout.write(restore_complex(pf, rec, ref_lookup, rg_tag) + "\n")
        n_restored += 1

    # Every layer must be exhausted. A leftover record means it pointed past the
    # end of the pBAM, i.e. the diff and the pBAM are not a matching pair --
    # almost always a per-region diff run against the wrong region's pBAM.
    rc = 0
    for r in readers:
        if r.line is not None:
            sys.stderr.write("ERROR: %s has records past the end of the pBAM "
                             "(next at line %d, pBAM has %d). The diff and pBAM "
                             "are not a matching pair.\n" % (r.path, r.line, lineno))
            rc = 1
        r.close()

    sys.stderr.write(
        "done.\n"
        "  pBAM lines:   %d\n"
        "  restored:     %d  (S=%d C=%d U=%d)\n"
        "  passed through: %d  (reference-only)\n"
        % (lineno, n_restored, n_cls["S"], n_cls["C"], n_cls["U"], n_through))
    if lineno != n_restored + n_through:
        sys.stderr.write("  WARNING: lines != restored + passed. Every pBAM line "
                         "must be emitted exactly once.\n")
        rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
