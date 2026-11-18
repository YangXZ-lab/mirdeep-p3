#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""Mature-boundary audit and correction.

The annotation step takes the mature interval from the prediction and then
adjusts it per arm (see ``utils/coords.py``).  That adjustment can leave the
interval one base out of frame at its 5' end relative to the reads that
actually support the locus.  The table still looks plausible, so nothing
complains.

This module measures the discrepancy from the bowtie alignment and re-anchors
the interval onto the read-supported position when the evidence justifies it.

Interval convention
-------------------
``basic-info`` carries each feature as ``start<TAB>end`` with a 1-based start
and an exclusive end, i.e. ``end - start`` equals the sequence length.  This is
not a guess: ``Nip-MIR166a`` carries ``3487876``/``3487897`` for a 21 nt mature.
Everything below is 0-based half-open internally, which matches that convention
directly.

The precursor sequence is written in transcript (mature) orientation: co-linear
with the genomic interval on ``+``, reverse complemented on ``-``.

Read offsets
------------
The bowtie reference is ``1 nt upstream flank + mature + 1 nt downstream
flank`` in transcript orientation.  ``offset == 1`` therefore means the read
starts on mature base 1, and a dominant offset of ``k`` asks for a shift of
``k - 1`` bases.

Correction is a shift, not a trim
---------------------------------
Moving only the 5' end is wrong as often as it is right.  ``Nip-MIR166a`` is one
base out of frame with the same length, so the interval has to be *shifted*;
``Gmo-MIR159a`` is one base too long, so it has to be *trimmed*.  The two cases
are indistinguishable from the offset alone, so both candidates are built and
the one that reproduces the dominant read sequence is accepted.  If neither
does, the locus is left alone.
"""

from typing import Dict, List, Optional, Tuple

# Correcting more than one base is not supported: at that distance the dominant
# offset is far more likely to be 5'-degradation than a boundary slip.
MAX_SHIFT = 1

# A corrected mature outside this range is treated as an artefact rather than a
# real correction.
MIN_MATURE_LEN = 18
MAX_MATURE_LEN = 26

OK = 'ok'
LOW_COVERAGE = 'low_coverage'
AMBIGUOUS = 'ambiguous'
SHIFT = 'shift'
SHIFT_UPSTREAM = 'shift_upstream'
SHIFT_BEYOND_LIMIT = 'shift_beyond_limit'
OUT_OF_RANGE = 'out_of_range'

# How the interval followed the reads.
ACTION_TRIM = 'trim'      # move the 5' end only
ACTION_SHIFT = 'shift'    # move the whole interval


def mature_rel_bounds(pre_start: int, pre_end: int, mature_start: int,
                      mature_end: int, strand: str) -> Tuple[int, int]:
    """
    Convert a genomic half-open mature interval to precursor-relative 0-based
    half-open bounds, in transcript orientation.

    ``pre_start`` / ``pre_end`` / ``mature_*`` follow the basic-info convention:
    1-based start, exclusive end.
    """
    if strand == '+':
        return mature_start - pre_start, mature_end - pre_start
    return pre_end - mature_end, pre_end - mature_start


def rel_to_genomic(pre_start: int, pre_end: int, rel_beg: int, rel_end: int,
                   strand: str) -> Tuple[int, int]:
    """Convert precursor-relative 0-based half-open bounds back to basic-info
    coordinates (1-based start, exclusive end)."""
    if strand == '+':
        return pre_start + rel_beg, pre_start + rel_end
    return pre_end - rel_end, pre_end - rel_beg


def offset_shift(offset: int) -> int:
    """
    Bases to move the mature 5' end downstream.

    ``offset == 1`` is canonical and asks for no movement.  A dominant offset of
    ``k`` means the reads start at mature-relative position ``k - 1``.
    """
    return offset - 1


def pick_dominant_offset(weights: Dict[int, float], min_weight: float,
                         min_fraction: float
                         ) -> Tuple[Optional[int], float, float]:
    """
    Pick the dominant 5' offset.

    Returns ``(offset, fraction, total_weight)``; ``offset`` is None when the
    locus has too little support to be judged.
    """
    total = sum(weights.values())
    if total < min_weight or total <= 0:
        return None, 0.0, total
    # Ties break towards the smaller offset, i.e. towards canonical, so a
    # perfectly bimodal locus is never corrected on a coin flip.
    offset, weight = max(weights.items(), key=lambda kv: (kv[1], -kv[0]))
    return offset, weight / total, total


def at_arm_boundary(rel_beg: int, rel_end: int, pre_len: int,
                    tolerance: int = 2) -> str:
    """
    Where on the hairpin the mature sits.

    ``'5p'`` when it starts at the precursor 5' end, ``'3p'`` when it ends at
    the precursor 3' end, ``'internal'`` otherwise.  Reported because a 3'-arm
    mature has a free 5' boundary and is therefore the likelier place for an
    annotation slip.
    """
    if rel_beg <= tolerance:
        return '5p'
    if pre_len - rel_end <= tolerance:
        return '3p'
    return 'internal'


def classify(*, weights: Dict[int, float], min_weight: float,
             min_fraction: float, max_shift: int = MAX_SHIFT,
             allow_extend: bool = False
             ) -> Dict[str, object]:
    """
    Decide whether a locus is a candidate for correction.

    ``status`` is one of the module constants; ``shift`` is the requested
    movement in bases and is only meaningful for status ``SHIFT``.  A ``SHIFT``
    verdict means "the reads ask for this shift"; whether a valid interval
    exists is decided by :func:`resolve_candidate`.

    A shift upstream (offset 0, reads starting before the annotated mature) is
    reported as ``SHIFT_UPSTREAM`` and left alone unless ``allow_extend`` is
    set.  Extending the 5' end means claiming a base the annotation deliberately
    excluded, and the 1 nt flank exists precisely to absorb such reads; treating
    them as a boundary call is a much stronger claim than trimming one base the
    reads do not support.
    """
    offset, fraction, total = pick_dominant_offset(weights, min_weight,
                                                   min_fraction)
    base = {'offset': offset, 'fraction': fraction, 'total': total, 'shift': 0}
    if offset is None:
        base['status'] = LOW_COVERAGE
        return base
    if offset == 1:
        base['status'] = OK
        return base
    if fraction < min_fraction:
        base['status'] = AMBIGUOUS
        return base
    shift = offset_shift(offset)
    base['shift'] = shift
    if abs(shift) > max_shift:
        base['status'] = SHIFT_BEYOND_LIMIT
        return base
    if shift < 0 and not allow_extend:
        base['status'] = SHIFT_UPSTREAM
        return base
    base['status'] = SHIFT
    return base


def candidate_intervals(rel_beg: int, rel_end: int, shift: int
                        ) -> Dict[str, Tuple[int, int]]:
    """
    The two ways an interval can follow the reads.

    ``trim`` moves the 5' end only, which changes the length; ``shift`` moves
    the whole interval, which preserves it.  Which one is right cannot be read
    off the offset, so both are offered and the caller picks the one that
    reproduces the reads.
    """
    return {
        ACTION_TRIM: (rel_beg + shift, rel_end),
        ACTION_SHIFT: (rel_beg + shift, rel_end + shift),
    }


def resolve_candidate(*, pre_seq: str, rel_beg: int, rel_end: int, shift: int,
                      read_seq: str
                      ) -> Optional[Tuple[str, int, int, str]]:
    """
    Pick the candidate interval that reproduces ``read_seq``.

    Returns ``(action, rel_beg, rel_end, seq)`` or None when neither candidate
    matches the reads, leaves the precursor, or yields an implausible length.
    Requiring an exact match is what keeps this from inventing an interval that
    merely happens to be allowed.
    """
    pre_seq = pre_seq.upper().replace('U', 'T')
    read_seq = read_seq.upper().replace('U', 'T')
    pre_len = len(pre_seq)
    if not read_seq or pre_len == 0:
        return None

    matches: List[Tuple[str, int, int, str]] = []
    for action, (beg, end) in candidate_intervals(rel_beg, rel_end,
                                                  shift).items():
        if beg < 0 or end > pre_len or beg >= end:
            continue
        length = end - beg
        if length < MIN_MATURE_LEN or length > MAX_MATURE_LEN:
            continue
        seq = pre_seq[beg:end]
        if seq == read_seq:
            matches.append((action, beg, end, seq))

    if len(matches) != 1:
        # Zero matches: the reads disagree with both candidate framings.
        # Two matches: the two actions coincide, which cannot happen for a
        # non-zero shift, so treat it as unexpected rather than pick one.
        return None
    return matches[0]


def corrected_genomic(*, pre_start: int, pre_end: int, strand: str,
                      rel_beg: int, rel_end: int) -> Optional[Tuple[int, int]]:
    """Map a resolved interval back to basic-info coordinates."""
    g_start, g_end = rel_to_genomic(pre_start, pre_end, rel_beg, rel_end,
                                    strand)
    if g_start <= 0 or g_end <= 0 or g_start >= g_end:
        return None
    return g_start, g_end
