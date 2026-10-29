#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""Coordinate arithmetic shared by the prediction and annotation steps.

Both conversions below are easy to get subtly wrong, and both fail *silently*
when wrong: a misplaced star still produces a plausible-looking basic-info
table, and a mirrored mature interval still produces a plausible-looking BED.
They therefore live in one place and are covered by tests/test_coordinates.py.

Two interval conventions are used on purpose -- do not "unify" them, they
mirror the input formats they consume:

* compute_star works on precursor-relative half-open intervals [start, end),
  length == end - start, matching the ranges read from the .struc file.
* map_prediction works on genomic 1-based inclusive intervals [start, end],
  length == end - start + 1, matching what basic-info / BED carry.
"""

from typing import Dict, Optional, Tuple

STAR_LEN = 21


def compute_star(pre_seq: str, mature_seq: str, pre_start: int, pre_end: int,
                 strand: str, star_len: int = STAR_LEN
                 ) -> Optional[Tuple[str, int, int]]:
    """
    Locate the miRNA* (passenger strand) on the precursor.

    ``pre_seq`` is written in transcript (mature) orientation: it is the strand
    on which ``mature_seq`` is found.  It is therefore co-linear with the
    genomic coordinates only for the '+' strand.  For the '-' strand it is the
    reverse complement of that interval, so its 5' end maps to ``pre_end`` and
    the genomic axis has to be mirrored.

    Coordinates follow the same convention as the input ranges, i.e. the
    interval covers [start, end) and its length equals ``end - start``.

    Returns:
        ``(star_seq, star_beg, star_end)``, or ``None`` if ``mature_seq``
        cannot be placed at either end of the precursor.
    """
    pre_len = len(pre_seq)
    mature_len = len(mature_seq)
    ind = pre_seq.find(mature_seq)
    if ind < 0:
        return None

    at5 = abs(ind) <= 2                              # mature at the precursor 5' end
    at3 = abs(ind - (pre_len - mature_len)) <= 2     # mature at the precursor 3' end
    if not (at5 or at3):
        return None

    # miRNA* occupies the opposite arm of the hairpin
    star_seq = pre_seq[max(0, pre_len - star_len):] if at5 else pre_seq[:star_len]
    n = len(star_seq)

    # '+' : pre_seq 5' end == pre_start -> star on the precursor 3' arm sits HIGH
    # '-' : pre_seq 5' end == pre_end   -> star on the precursor 3' arm sits LOW
    star_at_high = at5 if strand == '+' else (not at5)
    if star_at_high:
        return star_seq.upper(), pre_end - n, pre_end
    return star_seq.upper(), pre_start, pre_start + n


def map_prediction(*, mature_beg: int, mature_end: int, star_beg: int,
                   star_end: int, mature_arm: str, strand: str,
                   excise_beg: int, excise_end: int) -> Dict[str, int]:
    """
    Convert miRDeep-style precursor-relative coordinates to genomic ones.

    The relative coordinates are 1-based and follow transcript (mature)
    orientation; ``excise_beg`` / ``excise_end`` are the genomic bounds of the
    excised precursor window, also 1-based inclusive.

    The returned intervals are 1-based inclusive.  Two arm-dependent quirks are
    preserved from the original implementation because downstream BED files are
    calibrated against them:

    * ``mature_arm == 'first'`` shifts the mature 5' end by -1 and then extends
      the 3' end by +1, which compensates the 2-nt 3' overhang of the duplex.
    * ``mature_arm == 'second'`` shifts the mature 5' end by -1 on '+' and
      extends the 3' end by +1 on '-'.

    Raises:
        ValueError: with a short reason (used verbatim in the caller's warning)
            when a coordinate lands at or below zero.  The caller decides
            whether to drop the prediction.
    """
    if strand == '+':
        mature_bg_abs = excise_beg + mature_beg - 1
        mature_end_abs = excise_beg + mature_end - 1
        star_bg_abs = excise_beg + star_beg - 1
        star_end_abs = excise_beg + star_end - 1

        # Adjust precursor boundaries based on mature_arm
        if mature_arm == 'first':
            mature_bg_abs -= 1
            mature_end_abs -= 1
            pre_bg_abs = mature_bg_abs
            pre_end_abs = star_end_abs
        elif mature_arm == 'second':
            pre_bg_abs = star_bg_abs - 1
            pre_end_abs = mature_end_abs
        else:
            pre_bg_abs = excise_beg
            pre_end_abs = excise_end

        # Final adjustments
        if mature_arm == 'first':
            mature_end_abs += 1
        elif mature_arm == 'second':
            mature_bg_abs -= 1

    else:
        mature_bg_abs = excise_end - mature_end
        mature_end_abs = excise_end - mature_beg
        star_bg_abs = excise_end - star_end
        star_end_abs = excise_end - star_beg

        # Adjust precursor boundaries based on mature_arm
        if mature_arm == 'first':
            mature_bg_abs += 1
            mature_end_abs += 1
            pre_bg_abs = star_end_abs - 1
            pre_end_abs = mature_end_abs
        elif mature_arm == 'second':
            pre_bg_abs = mature_bg_abs
            pre_end_abs = star_end_abs
        else:
            pre_bg_abs = excise_beg
            pre_end_abs = excise_end

        # Final adjustments
        if mature_arm == 'first':
            mature_bg_abs -= 1
            pre_bg_abs = star_bg_abs
        elif mature_arm == 'second':
            mature_end_abs += 1
            pre_end_abs += 1

    # Verify coordinates are positive
    if (mature_bg_abs <= 0 or mature_end_abs <= 0
            or star_bg_abs <= 0 or star_end_abs <= 0
            or pre_bg_abs <= 0 or pre_end_abs <= 0):
        raise ValueError("negative or zero positions")

    if pre_bg_abs > pre_end_abs:
        pre_bg_abs, pre_end_abs = pre_end_abs, pre_bg_abs
    if mature_bg_abs > mature_end_abs:
        mature_bg_abs, mature_end_abs = mature_end_abs, mature_bg_abs
    if star_bg_abs > star_end_abs:
        star_bg_abs, star_end_abs = star_end_abs, star_bg_abs

    return {
        'mature_bg_abs': mature_bg_abs,
        'mature_end_abs': mature_end_abs,
        'star_bg_abs': star_bg_abs,
        'star_end_abs': star_end_abs,
        'pre_bg_abs': pre_bg_abs,
        'pre_end_abs': pre_end_abs,
    }
