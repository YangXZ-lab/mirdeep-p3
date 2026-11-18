#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""Unit tests for src/utils/mature_boundary.py.

Runs with the standard library only:

    python3 tests/test_mature_boundary.py

The fixture at the end is the real ``Nip-MIR166a`` record from a rice dataset:
its annotated mature is one base out of frame with the reads, which are the
miRBase mature.  It pins the interval convention (1-based start, exclusive end)
and the shift-versus-trim decision, both of which a wrong implementation gets
backwards without failing.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, "src"))

from utils.mature_boundary import (  # noqa: E402
    ACTION_SHIFT, ACTION_TRIM, AMBIGUOUS, LOW_COVERAGE, OK, SHIFT,
    SHIFT_BEYOND_LIMIT, SHIFT_UPSTREAM, at_arm_boundary, candidate_intervals, classify,
    corrected_genomic, mature_rel_bounds, offset_shift, pick_dominant_offset,
    rel_to_genomic, resolve_candidate,
)

FAILURES = []


def check(label, got, want):
    if got != want:
        FAILURES.append(f"{label}\n  got:  {got}\n  want: {want}")


# A 100-mer with no long internal repeat, so slice assertions stay meaningful.
PRE = ("ACGTTGCAAGTCCTAGGATCCGTAACGTTAGGCATCCGATTACAGGTCCAT"
       "GGATCCTTAACGGCTAAGCTTCGAGTCAGGCATTAACGGTTCAGGCATA")
check("fixture length", len(PRE), 100)

# Intervals follow basic-info: 1-based start, exclusive end.
PLUS_PRE_START, PLUS_PRE_END = 1000, 1100      # length 100
PLUS_MAT_START, PLUS_MAT_END = 1005, 1026      # length 21
MINUS_PRE_START, MINUS_PRE_END = 2000, 2100
MINUS_MAT_START, MINUS_MAT_END = 2050, 2071


def test_rel_bounds():
    check("rel_bounds '+'",
          mature_rel_bounds(PLUS_PRE_START, PLUS_PRE_END,
                            PLUS_MAT_START, PLUS_MAT_END, "+"), (5, 26))
    check("rel_bounds '-'",
          mature_rel_bounds(MINUS_PRE_START, MINUS_PRE_END,
                            MINUS_MAT_START, MINUS_MAT_END, "-"), (29, 50))
    check("rel_bounds length equals sequence length",
          (lambda b: b[1] - b[0])(mature_rel_bounds(
              PLUS_PRE_START, PLUS_PRE_END, PLUS_MAT_START, PLUS_MAT_END,
              "+")), PLUS_MAT_END - PLUS_MAT_START)


def test_rel_genomic_round_trip():
    for strand, ps, pe, ms, me in (
            ("+", PLUS_PRE_START, PLUS_PRE_END, PLUS_MAT_START, PLUS_MAT_END),
            ("-", MINUS_PRE_START, MINUS_PRE_END, MINUS_MAT_START, MINUS_MAT_END)):
        beg, end = mature_rel_bounds(ps, pe, ms, me, strand)
        check(f"round trip '{strand}'",
              rel_to_genomic(ps, pe, beg, end, strand), (ms, me))


def test_offset_shift():
    check("shift offset 1", offset_shift(1), 0)
    check("shift offset 2", offset_shift(2), 1)
    check("shift offset 0", offset_shift(0), -1)
    check("shift offset 3", offset_shift(3), 2)


def test_pick_dominant():
    check("dominant too light",
          pick_dominant_offset({1: 4.0}, min_weight=10, min_fraction=0.6),
          (None, 0.0, 4.0))
    offset, frac, total = pick_dominant_offset({1: 30.0, 2: 70.0},
                                               min_weight=10, min_fraction=0.6)
    check("dominant offset", offset, 2)
    check("dominant fraction", round(frac, 3), 0.7)
    check("dominant total", total, 100.0)


def test_classify():
    check("classify canonical",
          classify(weights={1: 100.0}, min_weight=10,
                   min_fraction=0.6)["status"], OK)
    check("classify low coverage",
          classify(weights={2: 3.0}, min_weight=10,
                   min_fraction=0.6)["status"], LOW_COVERAGE)
    check("classify ambiguous",
          classify(weights={2: 55.0, 1: 45.0}, min_weight=10,
                   min_fraction=0.6)["status"], AMBIGUOUS)
    # A 50/50 split must not be corrected on a coin flip.
    check("classify tie breaks to canonical",
          classify(weights={1: 50.0, 2: 50.0}, min_weight=10,
                   min_fraction=0.6)["status"], OK)
    shift = classify(weights={2: 80.0, 1: 20.0}, min_weight=10,
                     min_fraction=0.6)
    check("classify shift status", shift["status"], SHIFT)
    check("classify shift amount", shift["shift"], 1)
    check("classify beyond limit",
          classify(weights={3: 90.0, 1: 10.0}, min_weight=10,
                   min_fraction=0.6)["status"], SHIFT_BEYOND_LIMIT)
    # Reads starting upstream of the mature must not silently extend it.
    check("classify upstream reported only",
          classify(weights={0: 80.0, 1: 20.0}, min_weight=10,
                   min_fraction=0.6)["status"], SHIFT_UPSTREAM)
    check("classify upstream when allowed",
          classify(weights={0: 80.0, 1: 20.0}, min_weight=10,
                   min_fraction=0.6, allow_extend=True)["status"], SHIFT)


def test_candidates():
    check("candidates at shift 1", candidate_intervals(5, 26, 1),
          {ACTION_TRIM: (6, 26), ACTION_SHIFT: (6, 27)})
    check("candidates at shift -1", candidate_intervals(5, 26, -1),
          {ACTION_TRIM: (4, 26), ACTION_SHIFT: (4, 25)})


def test_resolve_candidate():
    # A read that keeps the length means the whole interval moved.
    got = resolve_candidate(pre_seq=PRE, rel_beg=5, rel_end=26, shift=1,
                            read_seq=PRE[6:27])
    check("resolve shift", got, (ACTION_SHIFT, 6, 27, PRE[6:27]))
    # A read that is one base shorter means only the 5' end moved.
    got = resolve_candidate(pre_seq=PRE, rel_beg=5, rel_end=26, shift=1,
                            read_seq=PRE[6:26])
    check("resolve trim", got, (ACTION_TRIM, 6, 26, PRE[6:26]))
    # Reads that match neither candidate leave the locus alone.
    check("resolve no match",
          resolve_candidate(pre_seq=PRE, rel_beg=5, rel_end=26, shift=1,
                            read_seq="ACGTACGT"), None)
    # A candidate running off the precursor is rejected even if it matched.
    check("resolve off precursor",
          resolve_candidate(pre_seq=PRE, rel_beg=0, rel_end=20, shift=-1,
                            read_seq=PRE[0:20]), None)


def test_corrected_genomic():
    check("genomic '+'", corrected_genomic(pre_start=1000, pre_end=1100,
                                           strand="+", rel_beg=6, rel_end=27),
          (1006, 1027))
    # The transcript 5' end is the genomic HIGH coordinate on '-'.
    check("genomic '-'", corrected_genomic(pre_start=2000, pre_end=2100,
                                           strand="-", rel_beg=30, rel_end=50),
          (2050, 2070))


def test_arm_boundary():
    check("arm internal", at_arm_boundary(5, 26, 100), "internal")
    check("arm 5p", at_arm_boundary(1, 22, 100), "5p")
    check("arm 3p", at_arm_boundary(78, 100, 100), "3p")


# ----------------------------------------------------------------------
# Real record: Nip-MIR166a, rice.  The annotated mature keeps one extra base
# at the 5' end and is one base short at the 3' end relative to the reads.
# ----------------------------------------------------------------------
PRE166 = ("ATTTTTATGGCAGTTGAGGGGAATGTTGTCTGGTCCGAGACCTAACACCGGGCGGAATGGCGGA"
          "TTCAGCTGCAGCTAAGCAAGCTAGGTGGGGGGTTTCGGACCAGGCTTCATTCCCCCCAACTCA"
          "ACTGTTTAGTT")
MIR166_ANNOTATED = "TTCGGACCAGGCTTCATTCCC"     # 21 nt, as written by 3.1.8b
MIR166_READS = "TCGGACCAGGCTTCATTCCCC"         # 21 nt, miRBase osa-miR166a-3p
MIR166_PRE_START, MIR166_PRE_END = 3487779, 3487917
MIR166_MAT_START, MIR166_MAT_END = 3487876, 3487897


def test_mir166a_regression():
    check("MIR166a fixture length", len(PRE166), 138)
    check("MIR166a annotated interval matches its sequence",
          PRE166[MIR166_MAT_START - MIR166_PRE_START:
                 MIR166_MAT_END - MIR166_PRE_START], MIR166_ANNOTATED)

    beg, end = mature_rel_bounds(MIR166_PRE_START, MIR166_PRE_END,
                                 MIR166_MAT_START, MIR166_MAT_END, "+")
    check("MIR166a rel bounds", (beg, end), (97, 118))

    verdict = classify(weights={2: 976202.0, 1: 615471.0}, min_weight=10,
                       min_fraction=0.60)
    check("MIR166a is a shift candidate", verdict["status"], SHIFT)
    check("MIR166a shift", verdict["shift"], 1)

    resolved = resolve_candidate(pre_seq=PRE166, rel_beg=beg, rel_end=end,
                                 shift=verdict["shift"], read_seq=MIR166_READS)
    check("MIR166a resolves to the read",
          resolved, (ACTION_SHIFT, 98, 119, MIR166_READS))
    check("MIR166a corrected coordinates",
          corrected_genomic(pre_start=MIR166_PRE_START,
                            pre_end=MIR166_PRE_END, strand="+",
                            rel_beg=98, rel_end=119),
          (3487877, 3487898))


def main():
    for test in (test_rel_bounds, test_rel_genomic_round_trip,
                 test_offset_shift, test_pick_dominant, test_classify,
                 test_candidates, test_resolve_candidate,
                 test_corrected_genomic, test_arm_boundary,
                 test_mir166a_regression):
        test()
    if FAILURES:
        print(f"FAILED {len(FAILURES)} check(s):")
        for failure in FAILURES:
            print(" " + failure)
        return 1
    print("test_mature_boundary: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
