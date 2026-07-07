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
"""
complex_variants.py — reference-free extraction of a complex read's full edit set from
CIGAR + MD + SEQ.

Purpose: give partition_diff enough to classify a complex read against the allowed list
(every SNV and every indel must be allowed for LP), and give pbam2bam enough to rebuild
the read's SEQ from the reference (subs to overwrite, inserted/soft-clip bases to splice
in). Indels are emitted in VCF-anchored (POS, REF, ALT) form, but WITHOUT left-alignment
(there is no reference here), so an aligner-placed indel that bcftools would have shifted
will simply miss the allowed list and route to MP -- fail-closed and safe.

A read is flagged force_mp when it carries content that can never be a catalogued allele
or that we cannot anchor: soft clips, unknown CIGAR ops, or an indel with no anchor base
(leading insertion, or an indel immediately after an N gap).
"""


def parse_cigar(cigar):
    ops = []
    num = ""
    for c in cigar:
        if c.isdigit():
            num += c
        else:
            ops.append((int(num), c))
            num = ""
    return ops


def parse_md(md):
    """Tokenize an MD string into ('=', n) match-runs, ('X', refbase), ('^', refbases)."""
    if md.startswith("MD:Z:"):
        md = md[5:]
    toks = []
    i = 0
    while i < len(md):
        c = md[i]
        if c.isdigit():
            j = i
            while j < len(md) and md[j].isdigit():
                j += 1
            toks.append(("=", int(md[i:j])))
            i = j
        elif c == "^":
            j = i + 1
            while j < len(md) and md[j].isalpha():
                j += 1
            toks.append(("^", md[i + 1:j]))
            i = j
        elif c.isalpha():
            toks.append(("X", c))
            i += 1
        else:
            i += 1
    return toks


class MDConsumer:
    """Pull reference info one M-base or one deletion at a time, in CIGAR order."""
    def __init__(self, md):
        self.toks = parse_md(md)
        self.ti = 0
        self.run = 0  # matches remaining in the current '=' run

    def next_mbase(self):
        """Return ('match', None) or ('mismatch', refbase) for one M/=/X base."""
        if self.run > 0:
            self.run -= 1
            return ("match", None)
        while self.ti < len(self.toks):
            kind, val = self.toks[self.ti]
            if kind == "=":
                self.ti += 1
                if val > 0:
                    self.run = val - 1
                    return ("match", None)
                continue  # zero-length run
            if kind == "X":
                self.ti += 1
                return ("mismatch", val)
            raise ValueError("MD: hit deletion while reading an aligned base")
        raise ValueError("MD: exhausted while reading an aligned base")

    def next_deletion(self):
        """At a CIGAR D op: consume the matching ^ token, return the deleted ref bases."""
        while self.ti < len(self.toks):
            kind, val = self.toks[self.ti]
            if kind == "=" and val == 0:
                self.ti += 1
                continue
            if kind == "^":
                self.ti += 1
                return val
            raise ValueError("MD: expected a deletion token, got %r" % ((kind, val),))
        raise ValueError("MD: exhausted at a deletion")


def extract(pos, cigar, seq, md, qual=None, min_bq=0):
    """
    Walk CIGAR+MD+SEQ. Returns a dict:
      subs    : list of (gpos, ref, alt)             M-portion SNVs (BQ-filtered)
      indels  : list of (gpos, ref, alt)             VCF-anchored, as-placed (no left-align)
      inserts : list of (read_offset, bases)         inserted bases, for reconstruction
      clips   : list of (read_offset, bases)         soft-clip bases, for reconstruction
      force_mp: bool                                  must route to MP regardless of alleles
    gpos is 1-based genomic. read_offset is 0-based into the original SEQ.
    """
    ops = parse_cigar(cigar)
    mdc = MDConsumer(md)
    gpos = pos
    ridx = 0
    subs, indels, inserts, clips = [], [], [], []
    last_ref_base = None
    force_mp = False

    for (L, op) in ops:
        if op in "M=X":
            for _ in range(L):
                kind, refbase = mdc.next_mbase()
                readbase = seq[ridx]
                if kind == "match":
                    rb = readbase
                else:
                    rb = refbase
                    if min_bq <= 0 or qual is None or (ord(qual[ridx]) - 33) >= min_bq:
                        subs.append((gpos, rb, readbase))
                last_ref_base = rb
                gpos += 1
                ridx += 1
        elif op == "I":
            ins = seq[ridx:ridx + L]
            if last_ref_base is None:
                force_mp = True            # leading insertion: no anchor base
            else:
                indels.append((gpos - 1, last_ref_base, last_ref_base + ins))
            inserts.append((ridx, ins))
            ridx += L
        elif op == "D":
            delbases = mdc.next_deletion()
            if last_ref_base is None:
                force_mp = True            # deletion with no anchor (e.g. after N)
            else:
                indels.append((gpos - 1, last_ref_base + delbases, last_ref_base))
            gpos += L
            last_ref_base = delbases[-1]
        elif op == "N":
            gpos += L
            last_ref_base = None           # cannot anchor across an intron
        elif op == "S":
            clips.append((ridx, seq[ridx:ridx + L]))
            ridx += L
            force_mp = True                # clipped sequence is not a catalogued allele
        elif op == "H":
            pass                           # not present in SEQ
        else:
            force_mp = True                # unknown op
    return {
        "subs": subs,
        "indels": indels,
        "inserts": inserts,
        "clips": clips,
        "force_mp": force_mp,
    }


def rebuild_seq(pos, cigar, ref_lookup, chrom, extracted):
    """
    Reconstruct the read SEQ from reference + extracted edits. Used by the validator now
    and (the same logic) by pbam2bam later. ref_lookup(chrom, start0, length) -> str.
    """
    ops = parse_cigar(cigar)
    subs = {g: a for (g, r, a) in extracted["subs"]}
    inserts = dict((o, b) for (o, b) in extracted["inserts"])
    clips = dict((o, b) for (o, b) in extracted["clips"])
    out = []
    gpos = pos
    ridx = 0
    for (L, op) in ops:
        if op in "M=X":
            ref = ref_lookup(chrom, gpos - 1, L)
            for k in range(L):
                g = gpos + k
                out.append(subs.get(g, ref[k]))
            gpos += L
            ridx += L
        elif op == "I":
            out.append(inserts.get(ridx, ""))
            ridx += L
        elif op == "D":
            gpos += L
        elif op == "N":
            gpos += L
        elif op == "S":
            out.append(clips.get(ridx, ""))
            ridx += L
        elif op == "H":
            pass
    return "".join(out)
