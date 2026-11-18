#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""Unit tests for src/utils/genome_fasta.py and src/utils/precursor_repair.py.

Standard library only:

    python3 tests/test_precursor_repair.py

The genome fetcher is exercised against a synthetic FASTA written to a temp
directory (multi-line, both strands), because the failure mode that matters --
an off-by-one in the byte arithmetic -- is invisible on a single-line record.
The repair logic is exercised with a fake fold, so the tests do not need
RNAfold installed; the real folding is covered by the end-to-end run.
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, "src"))

from utils.coords import compute_star  # noqa: E402
from utils.genome_fasta import GenomeFasta  # noqa: E402
from utils.precursor_repair import (  # noqa: E402
    CONTAINED, EXTENDED, MAX_EXTENSION, NO_FOLD, SEQ_MISMATCH,
    STAR_UNRESOLVED, TOO_LONG, added_bases, apply_repair, extension_needed,
    fold, repair_locus, reverse_complement, stem_window, strip_energy,
    transcript_offset,
)

FAILURES = []


def check(label, got, want):
    if got != want:
        FAILURES.append(f"{label}: got {got!r}, want {want!r}")


def check_true(label, condition):
    if not condition:
        FAILURES.append(f"{label}: expected true")


# --------------------------------------------------------------------------
# pure helpers
# --------------------------------------------------------------------------
check("revcomp", reverse_complement("ACGTN"), "NACGT")
check("revcomp roundtrip", reverse_complement(reverse_complement("ACGTTGCA")),
      "ACGTTGCA")

# mature already inside -> no work
check("no extension", extension_needed(100, 160, 120, 141), None)
check("flush both ends", extension_needed(100, 160, 100, 160), None)
# overhang at the high end (the Nip-MIR166a shape)
check("extend high", extension_needed(100, 160, 140, 161), (100, 161))
# overhang at the low end (the Nip-MIR166b shape)
check("extend low", extension_needed(100, 160, 99, 120), (99, 160))
# overhang on both sides at once
check("extend both", extension_needed(100, 160, 99, 161), (99, 161))

check("added high", added_bases(100, 160, 100, 161), (0, 1))
check("added low", added_bases(100, 160, 99, 160), (1, 0))
check("added both", added_bases(100, 160, 99, 161), (1, 1))

# where the old interval sits inside the widened one, in transcript orientation
check("offset plus", transcript_offset(100, 161, 100, 160, "+"), (0, 60))
check("offset plus low", transcript_offset(99, 160, 100, 160, "+"), (1, 61))
check("offset minus", transcript_offset(100, 161, 100, 160, "-"), (1, 61))
check("offset minus low", transcript_offset(99, 160, 100, 160, "-"), (0, 60))

check("stem window", stem_window(100, 160, 10_000), (80, 180))
check("stem window clamps start", stem_window(5, 60, 10_000), (1, 80))
check("stem window clamps end", stem_window(100, 160, 170), (80, 170))

check("strip energy", strip_energy("((((....)))) (-12.30)"), "((((....))))")
check("strip energy positive", strip_energy("((..)) ( 0.00)"), "((..))")
check("strip energy bare", strip_energy("((..))"), "((..))")


# --------------------------------------------------------------------------
# genome random access
# --------------------------------------------------------------------------
def write_fasta(path, records, width=7):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name} some description\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


CHR = "ACGTACGTACGTACGTACGT" + "TTGGCCAA" * 10 + "GGTTAACC" * 5   # 140 nt
OTHER = "AAAACCCCGGGGTTTT" * 4                                      # 64 nt
workdir = tempfile.mkdtemp(prefix="mirdp-precursor-")
genome_path = os.path.join(workdir, "mini.fa")
write_fasta(genome_path, [("chrX", CHR), ("chrY", OTHER)])

with GenomeFasta(genome_path) as g:
    check("records", sorted(g.records), ["chrX", "chrY"])
    check("length", g.length("chrX"), len(CHR))
    check("fetch full", g.fetch("chrX", 0, len(CHR)), CHR)
    check("fetch plus slice", g.fetch("chrX", 20, 30), CHR[20:30])
    check("fetch minus slice", g.fetch("chrX", 20, 30, "-"),
          reverse_complement(CHR[20:30]))
    check("fetch crosses a line", g.fetch("chrX", 5, 16), CHR[5:16])
    check("fetch last base", g.fetch("chrX", len(CHR) - 1, len(CHR)),
          CHR[-1:])
    check("fetch past end", g.fetch("chrX", 130, 141), None)
    check("fetch unknown chrom", g.fetch("chrZ", 0, 10), None)
    check("fetch empty range", g.fetch("chrX", 10, 10), None)

# a second open of the same file must reuse the cached index
with GenomeFasta(genome_path) as g2:
    check("cached index still correct", g2.fetch("chrY", 0, 8), OTHER[:8])


# --------------------------------------------------------------------------
# repair decisions
# --------------------------------------------------------------------------
def make_fake_genome(pre_start, pre_seq, mature_start, mature_end, strand,
                     flank=40):
    """Build a chromosome where `pre_seq` sits at ``pre_start``.

    The precursor is embedded with ``flank`` arbitrary bases on each side so a
    widening has something to read.
    """
    left = "ACGT" * (flank // 4 + 1)
    right = "TGCA" * (flank // 4 + 1)
    chrom = left[:flank] + pre_seq + right[:flank]
    return chrom


def fake_fold(seq):
    return "." * len(seq)


# A realistic shape: mature sits at the precursor 3' end (plus strand), exactly
# like Nip-MIR166a.  Then the mature is shifted one base downstream so it
# overhangs the precursor end by one.
PRE = ("GGTTTCGGACCAGGCTTCATTCCCCCCAACTCAACTGTTTAGTT"       # 43
       "ACGTACGTACGTACGTACGTA")                             # 21 -> 64 nt
PRE_START = 40
MAT_END = PRE_START + len(PRE)
MAT_START = MAT_END - 21
MATURE = PRE[-21:]
chrom = make_fake_genome(PRE_START, PRE, MAT_START, MAT_END, "+",
                         flank=40)


def fetch_plus(start, end, strand):
    if start < 0 or end > len(chrom) or start >= end:
        return None
    seq = chrom[start:end]
    return seq if strand == "+" else reverse_complement(seq)


# baseline: nothing to do
plan = repair_locus(pre_start=PRE_START, pre_end=MAT_END,
                    mature_start=MAT_START, mature_end=MAT_END,
                    mature_seq=MATURE, strand="+", old_pre_seq=PRE,
                    old_stem_seq=None, old_stem_struc=None,
                    fetch=fetch_plus, chrom_len=len(chrom),
                    compute_star=compute_star, fold_fn=fake_fold)
check("contained status", plan["status"], CONTAINED)
check("contained no change", plan["new_pre_seq"], None)

# the MIR166a shape: mature moved one base downstream, precursor is now 1 short
shift_start, shift_end = MAT_START + 1, MAT_END + 1
shifted_mature = chrom[shift_start:shift_end]
plan = repair_locus(pre_start=PRE_START, pre_end=MAT_END,
                    mature_start=shift_start, mature_end=shift_end,
                    mature_seq=shifted_mature, strand="+", old_pre_seq=PRE,
                    old_stem_seq=None, old_stem_struc=None,
                    fetch=fetch_plus, chrom_len=len(chrom),
                    compute_star=compute_star, fold_fn=fake_fold)
check("extended status", plan["status"], EXTENDED)
check("extended high by one", plan["added_high"], 1)
check("extended low not touched", plan["added_low"], 0)
check("extended length", len(plan["new_pre_seq"]), len(PRE) + 1)
check_true("old precursor still a prefix", plan["new_pre_seq"].startswith(PRE))
check_true("mature now inside", shifted_mature in plan["new_pre_seq"])
check("new pre end", plan["new_pre_end"], MAT_END + 1)
# the star sits on the opposite arm, so widening the 3' end must not move it
check("star unmoved start", plan["new_star_start"], PRE_START)
check("star length", plan["new_star_end"] - plan["new_star_start"], 21)

# widening beyond the limit is refused, not attempted
plan = repair_locus(pre_start=PRE_START, pre_end=MAT_END,
                    mature_start=PRE_START - (MAX_EXTENSION + 1),
                    mature_end=MAT_END, mature_seq=MATURE, strand="+",
                    old_pre_seq=PRE, old_stem_seq=None, old_stem_struc=None,
                    fetch=fetch_plus, chrom_len=len(chrom),
                    compute_star=compute_star, fold_fn=fake_fold)
check("too long", plan["status"], TOO_LONG)

# a genome that does not reproduce the stored precursor must not be edited
plan = repair_locus(pre_start=PRE_START, pre_end=MAT_END,
                    mature_start=shift_start, mature_end=shift_end,
                    mature_seq=shifted_mature, strand="+",
                    old_pre_seq=PRE[:-1] + "A" if PRE[-1] != "A" else PRE[:-1] + "C",
                    old_stem_seq=None, old_stem_struc=None,
                    fetch=fetch_plus, chrom_len=len(chrom),
                    compute_star=compute_star, fold_fn=fake_fold)
check("sequence mismatch", plan["status"], SEQ_MISMATCH)

# RNAfold failure leaves the row alone
plan = repair_locus(pre_start=PRE_START, pre_end=MAT_END,
                    mature_start=shift_start, mature_end=shift_end,
                    mature_seq=shifted_mature, strand="+", old_pre_seq=PRE,
                    old_stem_seq=None, old_stem_struc=None,
                    fetch=fetch_plus, chrom_len=len(chrom),
                    compute_star=compute_star, fold_fn=lambda s: None)
check("fold failure", plan["status"], NO_FOLD)

# minus strand, overhang at the low genomic end (the Nip-MIR166b shape).
# On '-' the low genomic coordinate is the transcript 3' end, so a mature that
# runs past pre_start is overhanging at its own 3' end -- exactly the case the
# rice MIR166b record shows.
pre_start_minus = PRE_START
pre_end_minus = PRE_START + len(PRE)
pre_seq_minus = reverse_complement(chrom[pre_start_minus:pre_end_minus])
# mature genomic [39, 60) -> 21 nt, one base below the precursor start
mat_start_minus, mat_end_minus = pre_start_minus - 1, pre_start_minus + 20
mature_minus = reverse_complement(chrom[mat_start_minus:mat_end_minus])
check("minus mature length", len(mature_minus), 21)
plan = repair_locus(pre_start=pre_start_minus, pre_end=pre_end_minus,
                    mature_start=mat_start_minus, mature_end=mat_end_minus,
                    mature_seq=mature_minus, strand="-",
                    old_pre_seq=pre_seq_minus,
                    old_stem_seq=None, old_stem_struc=None,
                    fetch=fetch_plus, chrom_len=len(chrom),
                    compute_star=compute_star, fold_fn=fake_fold)
check("minus status", plan["status"], EXTENDED)
check("minus extended low", plan["added_low"], 1)
check("minus extended high", plan["added_high"], 0)
# On '-' the added base lands at the end of the transcript, so the old
# precursor stays intact as a prefix (on '+' extending the high end would do
# the same; extending the low end prepends instead).
check_true("minus keeps old precursor as prefix",
           plan["new_pre_seq"].startswith(pre_seq_minus))
check_true("minus mature now inside", mature_minus in plan["new_pre_seq"])
check("minus star to the high side", plan["new_star_start"],
      pre_end_minus - 21)


# --------------------------------------------------------------------------
# writing the repair back into a row
# --------------------------------------------------------------------------
LAYOUT = {
    'stem_beg': 6, 'stem_end': 7, 'stem_seq': 8, 'stem_struc': 9,
    'pre_beg': 10, 'pre_end': 11, 'pre_seq': 12, 'pre_struc': 13,
    'star_beg': 19, 'star_end': 20, 'star_seq': 21,
}
row = [str(i) for i in range(22)]
row[15], row[16], row[17] = "100", "160", "MATURE_KEEP"
row[14], row[0], row[1] = "NAME_KEEP", "ID_KEEP", "ACC_KEEP"
plan = {
    'new_pre_start': 99, 'new_pre_end': 161,
    'new_pre_seq': "PRE", 'new_pre_struc': "PSTRUC",
    'new_stem_start': 79, 'new_stem_end': 181,
    'new_stem_seq': "STEM", 'new_stem_struc': "SSTRUC",
    'new_star_seq': "STAR", 'new_star_start': 100, 'new_star_end': 121,
}
before = list(row)
changed = apply_repair(row, plan, LAYOUT)
check("columns changed", changed, 11)
check("mature start untouched", row[15], "100")
check("mature end untouched", row[16], "160")
check("mature seq untouched", row[17], "MATURE_KEEP")
check("name untouched", row[0], "ID_KEEP")
check("mir_id untouched", row[1], "ACC_KEEP")
check("pre start written", row[10], "99")
check("pre seq written", row[12], "PRE")
check("star seq written", row[21], "STAR")
check("only allowed columns moved",
      [i for i in range(22) if row[i] != before[i]],
      [6, 7, 8, 9, 10, 11, 12, 13, 19, 20, 21])


# --------------------------------------------------------------------------
if FAILURES:
    print(f"test_precursor_repair: {len(FAILURES)} FAILURE(S)")
    for f in FAILURES:
        print("  -", f)
    sys.exit(1)
print("test_precursor_repair: all checks passed")
