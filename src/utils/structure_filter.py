#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""Reject loci whose two arms cannot both sit on the precursor.

A miRNA precursor has to hold a mature on one arm and a star on the other, with
a loop between them.  When the precursor is instead a window cut tight around
the mature, ``compute_star`` still "succeeds" -- it just takes ``star_len`` bases
from the opposite end -- and on a 40 nt precursor that lands on top of the mature:

    precursor 3742249..3742289   40 nt
    mature    3742249..3742269   20 nt
    star      3742268..3742289   21 nt      <- starts inside the mature
                                          20 + 21 = 41 > 40

Nothing downstream fails; the record simply carries two overlapping arms and a
dot-bracket that describes neither.  This module provides the missing check.

The measure is the *net gap* between the two intervals, which is the loop length:

    star after mature :  gap = star_start - mature_end
    star before mature:  gap = mature_start - star_end
    overlapping       :  gap < 0

A real hairpin loop is 15-100+ nt in plants; this module defaults to a
permissive 10 nt.  Note what the threshold actually tests: whether the precursor
has room for two arms and a loop, not whether the star is expressed.  Stars are
usually inferred rather than observed, so this is a structural sanity check.

Raising it is not free: on the rice dataset, gap < 10 nt rejects 14 of 446 loci
(3.1%), all of them with a precursor of 40-50 nt that is exactly
mature + star + loop.  Those loci are not repairable by widening the precursor
either -- every one of them has its mature in the middle of the +/-20 nt window
rather than on an arm, so ``compute_star`` on a wider precursor returns None.
"""

from typing import Dict, List, Optional, Tuple

#: Default minimum loop length between the two arms, in bases.
MIN_ARM_GAP = 10

#: Codes written to the report; also used as the ``status`` column.
KEPT = 'kept'
ARM_OVERLAP = 'arm_overlap'          # gap < 0, the arms cover each other
SHORT_ARM_GAP = 'short_arm_gap'      # 0 <= gap < min_gap
MISSING_PRECURSOR_STRUCTURE = 'missing_precursor_structure'
MISSING_STEM_STRUCTURE = 'missing_stem_structure'
MISSING_STAR = 'missing_star'
BAD_COORDINATES = 'bad_coordinates'

REJECTED = (ARM_OVERLAP, SHORT_ARM_GAP, MISSING_PRECURSOR_STRUCTURE,
            MISSING_STEM_STRUCTURE, MISSING_STAR, BAD_COORDINATES)

REPORT_HEADER = ['miRNA', 'status', 'reason', 'precursor_len', 'mature_len',
                 'star_len', 'arm_gap', 'precursor_start', 'precursor_end',
                 'mature_start', 'mature_end', 'star_start', 'star_end']

# 0-based column indices of the 22-column layout produced by annotation.
COL = {
    'name': 0, 'mir_id': 1, 'chrom': 4, 'strand': 5,
    'stem_start': 6, 'stem_end': 7, 'stem_seq': 8, 'stem_struc': 9,
    'pre_start': 10, 'pre_end': 11, 'pre_seq': 12, 'pre_struc': 13,
    'mat_name': 14, 'mat_start': 15, 'mat_end': 16, 'mat_seq': 17,
    'star_name': 18, 'star_start': 19, 'star_end': 20, 'star_seq': 21,
}


def arm_gap(mature_start: int, mature_end: int,
            star_start: int, star_end: int) -> int:
    """
    Net distance between the mature and star intervals.

    Positive is the loop length; zero means the two intervals abut; negative is
    the number of bases they overlap.  Strand-agnostic: it only looks at which
    interval lies to the low-coordinate side.
    """
    if star_start >= mature_end:
        return star_start - mature_end
    if mature_start >= star_end:
        return mature_start - star_end
    return -(min(mature_end, star_end) - max(mature_start, star_start))


def _structure_ok(seq: str, struc: str) -> bool:
    """A structure string is usable when it is non-empty and matches its sequence."""
    if not seq or not struc:
        return False
    return len(seq) == len(struc)


def assess(parts: List[str], min_gap: int = MIN_ARM_GAP) -> Dict[str, object]:
    """
    Decide whether one 22-column basic-info row survives the filter.

    Returns a dict with ``status``, ``reason`` (human readable, empty when kept)
    and the measurements used, so the report can be written without re-parsing.
    """
    out: Dict[str, object] = {
        'status': KEPT, 'reason': '',
        'precursor_len': 0, 'mature_len': 0, 'star_len': 0, 'arm_gap': '',
        'precursor_start': '', 'precursor_end': '',
        'mature_start': '', 'mature_end': '',
        'star_start': '', 'star_end': '',
    }

    def fail(status: str, reason: str) -> Dict[str, object]:
        out['status'] = status
        out['reason'] = reason
        return out

    try:
        pre_start = int(parts[COL['pre_start']])
        pre_end = int(parts[COL['pre_end']])
        mat_start = int(parts[COL['mat_start']])
        mat_end = int(parts[COL['mat_end']])
        star_start = int(parts[COL['star_start']])
        star_end = int(parts[COL['star_end']])
    except (ValueError, IndexError):
        return fail(BAD_COORDINATES, 'non-integer coordinates')

    out.update({
        'precursor_start': pre_start, 'precursor_end': pre_end,
        'mature_start': mat_start, 'mature_end': mat_end,
        'star_start': star_start, 'star_end': star_end,
        'precursor_len': pre_end - pre_start,
        'mature_len': mat_end - mat_start,
        'star_len': star_end - star_start,
    })

    if pre_end <= pre_start or mat_end <= mat_start or star_end <= star_start:
        return fail(BAD_COORDINATES, 'empty or inverted interval')

    # 1. both arms have to exist
    if not parts[COL['star_seq']].strip():
        return fail(MISSING_STAR, 'star sequence is empty')

    # 2. the precursor needs a structure, otherwise there is no hairpin at all
    if not _structure_ok(parts[COL['pre_seq']], parts[COL['pre_struc']]):
        return fail(MISSING_PRECURSOR_STRUCTURE,
                    'precursor structure missing or its length does not match '
                    'the precursor sequence')

    # 3. so does the stem window, which is what the structure plot draws
    if not _structure_ok(parts[COL['stem_seq']], parts[COL['stem_struc']]):
        return fail(MISSING_STEM_STRUCTURE,
                    'stem-window structure missing or its length does not match '
                    'the stem sequence')

    # 4. the two arms must not touch or overlap
    gap = arm_gap(mat_start, mat_end, star_start, star_end)
    out['arm_gap'] = gap
    if gap < 0:
        return fail(ARM_OVERLAP,
                    f'mature and star overlap by {-gap} nt '
                    f'({out["mature_len"]} + {out["star_len"]} > '
                    f'{out["precursor_len"]})')
    if gap < min_gap:
        return fail(SHORT_ARM_GAP,
                    f'loop between the arms is {gap} nt, below {min_gap}')

    return out


def filter_rows(rows: List[List[str]], min_gap: int = MIN_ARM_GAP
                ) -> Tuple[List[List[str]], List[Dict[str, object]]]:
    """
    Split rows into kept and rejected.

    Returns ``(kept_rows, report_rows)`` where ``report_rows`` is a dict per
    input row with its verdict, in input order.
    """
    kept: List[List[str]] = []
    report: List[Dict[str, object]] = []
    for parts in rows:
        verdict = assess(parts, min_gap=min_gap)
        entry = dict(verdict)
        entry['miRNA'] = parts[COL['name']] if parts else ''
        report.append(entry)
        if verdict['status'] == KEPT:
            kept.append(parts)
    return kept, report


def report_line(entry: Dict[str, object]) -> str:
    """Format one report entry as a TSV line (no trailing newline)."""
    return '\t'.join(str(entry.get(key, '')) for key in REPORT_HEADER)


def summarize(report: List[Dict[str, object]]) -> Dict[str, int]:
    """Count rows per status, in a stable order."""
    counts: Dict[str, int] = {}
    for entry in report:
        status = str(entry['status'])
        counts[status] = counts.get(status, 0) + 1
    return counts
