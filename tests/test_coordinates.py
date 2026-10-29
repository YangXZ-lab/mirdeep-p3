#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""Unit tests for src/utils/coords.py (O-05).

Runs with the standard library only:

    python3 tests/test_coordinates.py
    # or
    python3 -m pytest tests/test_coordinates.py

The four strand x arm combinations are checked against hand-computed values,
because those are exactly the combinations that a real run cannot distinguish
when the result is wrong: a mirrored star still looks like a miRNA.
"""

import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, "src"))

from utils.coords import compute_star, map_prediction  # noqa: E402

FAILURES = []


def check(label, got, want):
    if got != want:
        FAILURES.append(f"{label}\n  got:  {got}\n  want: {want}")


# ----------------------------------------------------------------------
# compute_star: precursor-relative, half-open intervals [start, end)
# ----------------------------------------------------------------------
# A pseudo-random 100-mer: a periodic sequence would let the 21-mer at one end
# also match near the other end, which defeats the point of these cases.
_RNG = random.Random(20261118)
PRE = "".join(_RNG.choice("ACGT") for _ in range(100))
PRE_LEN = len(PRE)
PRE_START, PRE_END = 1000, 1100        # half-open, length == end - start
MATURE = PRE[:21]                      # mature on the precursor 5' arm
STAR_HI = PRE[-21:]                    # star is the opposite (3') end
STAR_LO = PRE[:21]


def test_compute_star():
    # mature at the 5' end, '+' strand: star sits at the high coordinate.
    check("compute_star at5 '+'",
          compute_star(PRE, MATURE, PRE_START, PRE_END, "+"),
          (STAR_HI.upper(), 1079, 1100))
    # mature at the 5' end, '-' strand: pre_seq is the reverse complement of
    # the genomic interval, so its 5' end maps to pre_end and the star is low.
    check("compute_star at5 '-'",
          compute_star(PRE, MATURE, PRE_START, PRE_END, "-"),
          (STAR_HI.upper(), 1000, 1021))

    # mature at the 3' end of the precursor (mirror of the two cases above).
    mature3 = PRE[PRE_LEN - 21:]
    check("compute_star at3 '+'",
          compute_star(PRE, mature3, PRE_START, PRE_END, "+"),
          (STAR_LO.upper(), 1000, 1021))
    check("compute_star at3 '-'",
          compute_star(PRE, mature3, PRE_START, PRE_END, "-"),
          (STAR_LO.upper(), 1079, 1100))

    # Every reported interval must be as long as the sequence it describes.
    # 3.1.6c shipped a 21-nt star with a 20-nt interval (I-21 regression).
    for seq, beg, end in [
            compute_star(PRE, MATURE, PRE_START, PRE_END, "+"),
            compute_star(PRE, MATURE, PRE_START, PRE_END, "-"),
            compute_star(PRE, mature3, PRE_START, PRE_END, "+"),
            compute_star(PRE, mature3, PRE_START, PRE_END, "-")]:
        check("compute_star interval length", end - beg, len(seq))

    # Mature that sits in the middle of the precursor is not a miRNA candidate.
    check("compute_star middle", compute_star(PRE, PRE[30:51],
                                              PRE_START, PRE_END, "+"), None)
    # Mature that is not in the precursor at all.
    check("compute_star absent", compute_star(PRE, "TTTTTTTTTTTTTTTTTTTTT",
                                              PRE_START, PRE_END, "+"), None)


# ----------------------------------------------------------------------
# map_prediction: genomic, 1-based inclusive intervals [start, end]
# ----------------------------------------------------------------------
EXCISE_BEG, EXCISE_END = 1001, 1100    # 100-nt excised window, inclusive

# Relative coordinates of a 21-nt duplex: the arm called "first" has the mature
# at the low end, "second" at the high end (transcript orientation).
REL_FIRST = dict(mature_beg=5, mature_end=25, star_beg=40, star_end=60)
REL_SECOND = dict(mature_beg=40, mature_end=60, star_beg=5, star_end=25)

# Expected values, hand computed.  Mature is always one nt longer than its
# span because of the 3' overhang compensation in map_prediction.
EXPECTED = [
    ("'+' first arm", dict(strand="+", mature_arm="first", **REL_FIRST),
     dict(mature_bg_abs=1004, mature_end_abs=1025,
          star_bg_abs=1040, star_end_abs=1060,
          pre_bg_abs=1004, pre_end_abs=1060)),
    ("'+' second arm", dict(strand="+", mature_arm="second", **REL_SECOND),
     dict(mature_bg_abs=1039, mature_end_abs=1060,
          star_bg_abs=1005, star_end_abs=1025,
          pre_bg_abs=1004, pre_end_abs=1060)),
    ("'-' first arm", dict(strand="-", mature_arm="first", **REL_FIRST),
     dict(mature_bg_abs=1075, mature_end_abs=1096,
          star_bg_abs=1040, star_end_abs=1060,
          pre_bg_abs=1040, pre_end_abs=1096)),
    ("'-' second arm", dict(strand="-", mature_arm="second", **REL_SECOND),
     dict(mature_bg_abs=1040, mature_end_abs=1061,
          star_bg_abs=1075, star_end_abs=1095,
          pre_bg_abs=1040, pre_end_abs=1096)),
]


def test_map_prediction():
    for label, kwargs, want in EXPECTED:
        got = map_prediction(excise_beg=EXCISE_BEG, excise_end=EXCISE_END,
                             **kwargs)
        check(f"map_prediction {label}", got, want)

        # Structural invariants that must hold for every combination.
        check(f"{label}: mature span + 1 nt overhang",
              got["mature_end_abs"] - got["mature_bg_abs"] + 1,
              kwargs["mature_end"] - kwargs["mature_beg"] + 2)
        check(f"{label}: precursor contains mature",
              got["pre_bg_abs"] <= got["mature_bg_abs"]
              and got["pre_end_abs"] >= got["mature_end_abs"], True)
        # 3.1.6c regression: a star that overlaps the mature had the two arms
        # swapped (171/171 predictions in the test set were affected).
        check(f"{label}: star does not overlap mature",
              got["star_end_abs"] < got["mature_bg_abs"]
              or got["star_bg_abs"] > got["mature_end_abs"], True)
        check(f"{label}: all coordinates positive",
              min(got.values()) > 0, True)

    # '-' strand mirrors the genomic axis, so the arm order flips relative to
    # '+' : star below mature for the first arm, above it for the second.
    neg_first = map_prediction(excise_beg=EXCISE_BEG, excise_end=EXCISE_END,
                               strand="-", mature_arm="first", **REL_FIRST)
    check("'-' first arm: star below mature",
          neg_first["star_end_abs"] < neg_first["mature_bg_abs"], True)
    neg_second = map_prediction(excise_beg=EXCISE_BEG, excise_end=EXCISE_END,
                                strand="-", mature_arm="second", **REL_SECOND)
    check("'-' second arm: star above mature",
          neg_second["star_bg_abs"] > neg_second["mature_end_abs"], True)

    # Unknown arm: the precursor interval falls back to the excised window.
    fallback = map_prediction(excise_beg=EXCISE_BEG, excise_end=EXCISE_END,
                              strand="+", mature_arm="",
                              mature_beg=5, mature_end=25,
                              star_beg=40, star_end=60)
    check("unknown arm falls back to excised window",
          (fallback["pre_bg_abs"], fallback["pre_end_abs"]),
          (EXCISE_BEG, EXCISE_END))


def test_map_prediction_errors():
    # Excised window starting at 1 with a relative coordinate that lands on 0.
    try:
        map_prediction(excise_beg=1, excise_end=100, strand="+",
                       mature_arm="first", mature_beg=0, mature_end=20,
                       star_beg=30, star_end=50)
    except ValueError as exc:
        check("non-positive coordinate message", str(exc),
              "negative or zero positions")
    else:
        FAILURES.append("map_prediction did not raise on non-positive coords")


def test_interval_swap():
    # A precursor whose star/mature ranges are written in the opposite order
    # must still come back as start <= end.
    got = map_prediction(excise_beg=EXCISE_BEG, excise_end=EXCISE_END,
                         strand="+", mature_arm="second",
                         mature_beg=60, mature_end=40,
                         star_beg=25, star_end=5)
    for key in ("mature_bg_abs", "star_bg_abs", "pre_bg_abs"):
        end_key = key.replace("_bg_abs", "_end_abs")
        check(f"swap keeps {key} <= {end_key}",
              got[key] <= got[end_key], True)


def main():
    for test in (test_compute_star, test_map_prediction,
                 test_map_prediction_errors, test_interval_swap):
        test()
    if FAILURES:
        print(f"FAILED {len(FAILURES)} check(s):")
        for failure in FAILURES:
            print(" " + failure)
        return 1
    print("test_coordinates: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
