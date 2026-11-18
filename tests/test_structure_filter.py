#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""Unit tests for src/utils/structure_filter.py.

Standard library only:

    python3 tests/test_structure_filter.py

Runs on synthetic 22-column rows so each branch of the predicate is exercised
directly.  The numbers in the ``Nip-MIR12457a`` fixture are the real record --
a 40 nt precursor with a 20 nt mature and a 21 nt star, where 20 + 21 > 40 and
the star starts inside the mature.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, "src"))

from utils.structure_filter import (  # noqa: E402
    ARM_OVERLAP, BAD_COORDINATES, KEPT, MIN_ARM_GAP, MISSING_PRECURSOR_STRUCTURE,
    MISSING_STAR, MISSING_STEM_STRUCTURE, SHORT_ARM_GAP, arm_gap, assess,
    filter_rows, report_line, summarize,
)

FAILURES = []


def check(label, got, want):
    if got != want:
        FAILURES.append(f"{label}: got {got!r}, want {want!r}")


def check_true(label, condition):
    if not condition:
        FAILURES.append(f"{label}: expected true")


def row(*, name="MIRx", pre=(100, 200), mat=(100, 121), star=(179, 200),
        pre_struc=None, stem_struc=None, star_seq="ACGT"):
    """Build a 22-column row; defaults describe a normal hairpin."""
    parts = [""] * 22
    parts[0] = name
    parts[1] = "MIR_ID000000001"
    parts[4] = "chr01"
    parts[5] = "+"
    pre_len = pre[1] - pre[0]
    stem_len = pre_len + 40
    parts[6], parts[7] = str(pre[0] - 20), str(pre[1] + 20)
    parts[8] = "A" * stem_len
    parts[9] = stem_struc if stem_struc is not None else "." * stem_len
    parts[10], parts[11] = str(pre[0]), str(pre[1])
    parts[12] = "A" * pre_len
    parts[13] = pre_struc if pre_struc is not None else "." * pre_len
    parts[14] = f"{name.lower()}"
    parts[15], parts[16] = str(mat[0]), str(mat[1])
    parts[17] = "A" * (mat[1] - mat[0])
    parts[18] = f"{name.lower()}*"
    parts[19], parts[20] = str(star[0]), str(star[1])
    parts[21] = star_seq * max(1, (star[1] - star[0]) // 4)
    parts[21] = parts[21][:star[1] - star[0]]
    return parts


# --------------------------------------------------------------------------
# arm_gap: the four geometric cases
# --------------------------------------------------------------------------
check("gap star after", arm_gap(100, 121, 179, 200), 58)
check("gap star before", arm_gap(179, 200, 100, 121), 58)
check("gap abutting", arm_gap(100, 121, 121, 142), 0)
check("gap overlap by one", arm_gap(100, 121, 120, 141), -1)
check("gap overlap by five", arm_gap(100, 121, 116, 141), -5)

# the real Nip-MIR12457a geometry, minus strand: mature low, star high
check("gap MIR12457a", arm_gap(3742249, 3742269, 3742268, 3742289), -1)


# --------------------------------------------------------------------------
# assess: one branch at a time
# --------------------------------------------------------------------------
v = assess(row())
check("normal kept", v["status"], KEPT)
check("normal gap reported", v["arm_gap"], 58)
check("normal no reason", v["reason"], "")
check("normal precursor len", v["precursor_len"], 100)
check("normal mature len", v["mature_len"], 21)
check("normal star len", v["star_len"], 21)

# a loop right at the threshold survives, one base below does not
v = assess(row(pre=(100, 200), mat=(100, 130), star=(140, 200)))
check("gap exactly 10", v["arm_gap"], 10)
check("gap exactly 10 kept", v["status"], KEPT)
v = assess(row(pre=(100, 200), mat=(100, 130), star=(139, 200)))
check("gap 9", v["arm_gap"], 9)
check("gap 9 rejected", v["status"], SHORT_ARM_GAP)

# the threshold is configurable, and 0 keeps everything non-overlapping
v = assess(row(pre=(100, 200), mat=(100, 130), star=(139, 200)), min_gap=0)
check("gap 9 with min 0 kept", v["status"], KEPT)
v = assess(row(pre=(100, 200), mat=(100, 130), star=(139, 200)), min_gap=15)
check("gap 9 with min 15 rejected", v["status"], SHORT_ARM_GAP)

# the real record: 20 + 21 > 40
mir12457a = row(name="Nip-MIR12457a", pre=(3742249, 3742289),
                mat=(3742249, 3742269), star=(3742268, 3742289))
v = assess(mir12457a)
check("MIR12457a overlap", v["status"], ARM_OVERLAP)
check("MIR12457a gap", v["arm_gap"], -1)
check("MIR12457a precursor len", v["precursor_len"], 40)
check_true("MIR12457a reason mentions overlap", "overlap" in v["reason"])

# a precursor that cannot hold the two arms, seen through the gap instead
v = assess(row(pre=(0, 50), mat=(0, 21), star=(29, 50)))
check("50 nt precursor gap 8", v["arm_gap"], 8)
check("50 nt precursor rejected", v["status"], SHORT_ARM_GAP)

# missing pieces
v = assess(row(star_seq=""))
check("no star seq", v["status"], MISSING_STAR)

v = assess(row(pre=(100, 200), pre_struc=""))
check("no precursor structure", v["status"], MISSING_PRECURSOR_STRUCTURE)

v = assess(row(pre=(100, 200), pre_struc="." * 50))   # wrong length
check("precursor structure length mismatch", v["status"],
      MISSING_PRECURSOR_STRUCTURE)

v = assess(row(pre=(100, 200), stem_struc=""))
check("no stem structure", v["status"], MISSING_STEM_STRUCTURE)

v = assess(row(pre=(100, 200), stem_struc="." * 10))
check("stem structure length mismatch", v["status"], MISSING_STEM_STRUCTURE)

# an inverted interval is a coordinate bug, not a short loop
bad = row()
bad[16] = bad[15]
v = assess(bad)
check("empty mature", v["status"], BAD_COORDINATES)

bad = row()
bad[15] = "not-a-number"
v = assess(bad)
check("non-integer", v["status"], BAD_COORDINATES)


# --------------------------------------------------------------------------
# filter_rows / report formatting
# --------------------------------------------------------------------------
rows = [
    row(name="keep1"),
    mir12457a,
    row(name="short", pre=(0, 50), mat=(0, 21), star=(29, 50)),
    row(name="nostar", star_seq=""),
]
kept, report = filter_rows(rows)
check("kept count", len(kept), 1)
check("kept name", kept[0][0], "keep1")
check("report covers every input row", len(report), 4)
check("report order preserved", [e["miRNA"] for e in report],
      ["keep1", "Nip-MIR12457a", "short", "nostar"])

counts = summarize(report)
check("summary kept", counts.get(KEPT), 1)
check("summary overlap", counts.get(ARM_OVERLAP), 1)
check("summary short gap", counts.get(SHORT_ARM_GAP), 1)
check("summary missing star", counts.get(MISSING_STAR), 1)

line = report_line(report[1])
fields = line.split("\t")
check("report line column count", len(fields), 13)
check("report line name", fields[0], "Nip-MIR12457a")
check("report line status", fields[1], ARM_OVERLAP)
check("report line gap", fields[6], "-1")

# a kept row must survive untouched, byte for byte
kept, _ = filter_rows([row(name="untouched")])
check("kept row unchanged", kept[0], row(name="untouched"))

check("default threshold", MIN_ARM_GAP, 10)


# --------------------------------------------------------------------------
if FAILURES:
    print(f"test_structure_filter: {len(FAILURES)} FAILURE(S)")
    for f in FAILURES:
        print("  -", f)
    sys.exit(1)
print("test_structure_filter: all checks passed")
