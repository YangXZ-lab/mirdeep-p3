#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""Re-extend a precursor so it contains a corrected mature, and re-fold it.

A mature 5' correction can push the mature one base past the annotated
precursor end, because the hairpin was excised flush with the old (wrong) 3'
end.  Leaving it that way is internally inconsistent: the mature is no longer a
substring of the precursor, and anything that maps mature coordinates onto the
precursor sequence (the dot-bracket plot, the stored structure string) is now
mislabelled by one base.

This module closes that gap.  For a locus whose mature no longer fits, it

1. widens the precursor on the side(s) that are short, reading the extra bases
   from the genome;
2. verifies the widening only *added* bases -- the old precursor must still be
   present at the expected offset, otherwise the assembly disagrees and the
   locus is left alone;
3. re-folds the new precursor with RNAfold, and the extended +/-20 nt window
   as well (that window is what ``RNA_plot.py`` draws);
4. recomputes the star from the new precursor and reports whether it moved.

Only the columns that describe those two windows and the star are rewritten.
The mature interval, the names and the MIR_ID are never touched here.

Coordinate convention: every interval in this module is 0-based and half-open
(``[start, end)``, length == end - start), matching ``*-basic-info``.
"""

import subprocess
from typing import Callable, Dict, List, Optional, Tuple

#: A mature that needs more than this many extra bases on one side is not a
#: boundary problem any more; something else is wrong and we refuse to guess.
MAX_EXTENSION = 5

#: Width of the stem window stored alongside the precursor (``Precursor_20nt``).
STEM_FLANK = 20

CONTAINED = 'contained'          # mature already inside precursor, nothing to do
EXTENDED = 'extended'            # precursor widened, structures re-folded
TOO_LONG = 'extension_too_long'  # needs more than MAX_EXTENSION bases
OUT_OF_GENOME = 'out_of_genome'  # widening would leave the chromosome
SEQ_MISMATCH = 'sequence_mismatch'   # genome does not reproduce the precursor
STAR_UNRESOLVED = 'star_unresolved'  # widened, but compute_star found no arm
NO_FOLD = 'fold_failed'          # RNAfold did not return a structure
MATURE_MISSING = 'mature_missing'    # still not inside after widening

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_ENERGY_SUFFIX = __import__("re").compile(r"\s*\(\s*[-+]?[\d.]+\s*\)\s*$")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def extension_needed(pre_start: int, pre_end: int, mat_start: int,
                     mat_end: int) -> Optional[Tuple[int, int]]:
    """Return the smallest precursor interval containing the mature, or None."""
    if mat_start >= pre_start and mat_end <= pre_end:
        return None
    return min(pre_start, mat_start), max(pre_end, mat_end)


def added_bases(pre_start: int, pre_end: int, new_start: int,
                new_end: int) -> Tuple[int, int]:
    """(bases added at the low end, bases added at the high end)."""
    return pre_start - new_start, new_end - pre_end


def transcript_offset(new_start: int, new_end: int, old_start: int,
                      old_end: int, strand: str) -> Tuple[int, int]:
    """
    Where an old genomic interval sits inside a wider one, in the wider one's
    transcript-oriented sequence (0-based half-open).
    """
    if strand == '+':
        return old_start - new_start, old_end - new_start
    return new_end - old_end, new_end - old_start


def stem_window(pre_start: int, pre_end: int, chrom_len: int
                ) -> Tuple[int, int]:
    """The precursor +/-20 nt window, clamped exactly as the pipeline builds it."""
    return max(1, pre_start - STEM_FLANK), min(chrom_len, pre_end + STEM_FLANK)


def strip_energy(line: str) -> str:
    """Drop RNAfold's trailing ``( -12.30)`` from a dot-bracket line."""
    return _ENERGY_SUFFIX.sub("", line).strip()


def fold(seq: str, rnafold: str = "RNAfold") -> Optional[str]:
    """Fold one sequence and return its dot-bracket string, or None."""
    try:
        proc = subprocess.run([rnafold, "--noPS"], input=f">q\n{seq}\n",
                              capture_output=True, text=True, check=True)
    except (OSError, subprocess.CalledProcessError):
        return None
    lines = [l for l in proc.stdout.splitlines() if l.strip()]
    if len(lines) < 3:
        return None
    structure = strip_energy(lines[2])
    if len(structure) != len(seq):
        return None
    return structure


def repair_locus(*, pre_start: int, pre_end: int, mature_start: int,
                 mature_end: int, mature_seq: str, strand: str,
                 old_pre_seq: str, old_stem_seq: Optional[str],
                 old_stem_struc: Optional[str],
                 fetch: Callable[[int, int, str], Optional[str]],
                 chrom_len: int,
                 compute_star: Callable[..., Optional[Tuple[str, int, int]]],
                 fold_fn: Callable[[str], Optional[str]] = fold,
                 max_extension: int = MAX_EXTENSION) -> Dict[str, object]:
    """
    Plan (and fold) the repair for one locus.

    ``fetch`` is called as ``fetch(start, end, strand)`` and must return the
    transcript-oriented sequence.  ``compute_star`` has the signature of
    ``utils.coords.compute_star``; it is injected so this module stays free of
    import cycles and easy to test.

    Returns a dict with ``status`` plus the fields the caller may need; keys are
    always present so the caller can write a uniform report row.
    """
    out: Dict[str, object] = {
        'status': CONTAINED,
        'new_pre_start': pre_start, 'new_pre_end': pre_end,
        'new_pre_seq': None, 'new_pre_struc': None,
        'new_stem_start': None, 'new_stem_end': None,
        'new_stem_seq': None, 'new_stem_struc': None,
        'new_star_seq': None, 'new_star_start': None, 'new_star_end': None,
        'star_moved': None, 'added_low': 0, 'added_high': 0, 'note': '',
    }

    target = extension_needed(pre_start, pre_end, mature_start, mature_end)
    if target is None:
        return out
    new_start, new_end = target
    low, high = added_bases(pre_start, pre_end, new_start, new_end)
    out['added_low'], out['added_high'] = low, high
    out['new_pre_start'], out['new_pre_end'] = new_start, new_end

    if max(low, high) > max_extension:
        out['status'] = TOO_LONG
        out['note'] = f"needs {max(low, high)} nt, limit {max_extension}"
        return out
    if new_start < 0 or new_end > chrom_len:
        out['status'] = OUT_OF_GENOME
        out['note'] = f"would leave 0..{chrom_len}"
        return out

    new_pre_seq = fetch(new_start, new_end, strand)
    if new_pre_seq is None:
        out['status'] = OUT_OF_GENOME
        out['note'] = 'genome fetch returned nothing'
        return out

    # The widening must only have added bases: the old precursor has to still be
    # there at the offset the coordinates imply.  If it is not, the genome is a
    # different assembly and we must not touch the row.
    lo, hi = transcript_offset(new_start, new_end, pre_start, pre_end, strand)
    if new_pre_seq[lo:hi] != old_pre_seq.upper():
        out['status'] = SEQ_MISMATCH
        out['note'] = 'genome does not reproduce the stored precursor'
        return out
    if mature_seq.upper() not in new_pre_seq:
        out['status'] = MATURE_MISSING
        out['note'] = 'mature is not a substring of the widened precursor'
        return out

    new_pre_struc = fold_fn(new_pre_seq)
    if new_pre_struc is None:
        out['status'] = NO_FOLD
        out['note'] = 'RNAfold failed on the widened precursor'
        return out
    out['new_pre_seq'] = new_pre_seq
    out['new_pre_struc'] = new_pre_struc

    # The +/-20 nt window is derived from the precursor, so it moves with it.
    stem_start, stem_end = stem_window(new_start, new_end, chrom_len)
    stem_seq = fetch(stem_start, stem_end, strand)
    if stem_seq is not None:
        stem_struc = fold_fn(stem_seq)
        if stem_struc is not None and len(stem_struc) == len(stem_seq):
            out['new_stem_start'] = stem_start
            out['new_stem_end'] = stem_end
            out['new_stem_seq'] = stem_seq
            out['new_stem_struc'] = stem_struc
        else:
            out['note'] = 'stem window re-fold failed; kept the old window'

    star = compute_star(new_pre_seq, mature_seq.upper(), new_start, new_end,
                        strand)
    if star is None:
        out['status'] = STAR_UNRESOLVED
        out['note'] = (out['note'] + '; ' if out['note'] else '') + \
            'mature is not at either precursor arm, star not recomputed'
        return out
    star_seq, star_beg, star_end = star
    out['new_star_seq'] = star_seq
    out['new_star_start'] = star_beg
    out['new_star_end'] = star_end
    out['status'] = EXTENDED
    return out


def apply_repair(parts: List[str], plan: Dict[str, object],
                 layout: Dict[str, int]) -> int:
    """
    Write a repair plan into a 22-column basic-info row, in place.

    Returns the number of columns that actually changed.  Only the stem window,
    the precursor and the star are eligible; everything else is left alone.
    """
    changed = 0
    pairs = [
        (layout['stem_beg'], plan['new_stem_start']),
        (layout['stem_end'], plan['new_stem_end']),
        (layout['stem_seq'], plan['new_stem_seq']),
        (layout['stem_struc'], plan['new_stem_struc']),
        (layout['pre_beg'], plan['new_pre_start']),
        (layout['pre_end'], plan['new_pre_end']),
        (layout['pre_seq'], plan['new_pre_seq']),
        (layout['pre_struc'], plan['new_pre_struc']),
    ]
    for idx, value in pairs:
        if value is None or idx is None:
            continue
        text = str(value)
        if parts[idx] != text:
            parts[idx] = text
            changed += 1

    if plan.get('new_star_seq') is not None:
        star_pairs = [
            (layout.get('star_seq'), plan['new_star_seq']),
            (layout.get('star_beg'), plan['new_star_start']),
            (layout.get('star_end'), plan['new_star_end']),
        ]
        for idx, value in star_pairs:
            if idx is None or value is None:
                continue
            text = str(value)
            if parts[idx] != text:
                parts[idx] = text
                changed += 1
    return changed
