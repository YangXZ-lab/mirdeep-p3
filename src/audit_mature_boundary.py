#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""Audit the mature 5' boundary of a basic-info table against its alignment.

Reads per-locus 5' offsets from ``*.mature.aln`` (bowtie, reference =
1 nt flank + mature + 1 nt flank), keeps only reads that reach the 3' end, and
reports every locus whose dominant offset is not canonical.

With ``--apply`` the table is rewritten with the corrected mature interval.
A correction is only accepted when it reproduces the dominant read sequence
inside the precursor, so the output interval is the one the reads support rather
than a plausible-looking neighbour.

Inputs:
  -i:   basic-info table (22 or 24 columns)
  -aln: alignment files (comma separated) or -aln-dir (globbed *.mature.aln)
  -o:   audit report (TSV)
  -oc:  corrected basic-info (only written with --apply)
"""

import argparse
import glob
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils.mature_boundary import (  # noqa: E402
    AMBIGUOUS, LOW_COVERAGE, OK, OUT_OF_RANGE, SHIFT, SHIFT_BEYOND_LIMIT,
    SHIFT_UPSTREAM, MAX_SHIFT, at_arm_boundary, classify, corrected_genomic,
    mature_rel_bounds, resolve_candidate,
)

REPORT_HEADER = ['miRNA', 'status', 'dominant_offset', 'dominance', 'weight',
                 'arm', 'action', 'old_start', 'old_end', 'old_seq',
                 'new_start', 'new_end', 'new_seq']

# Column layouts, 0-based.  Two windows matter and they are not the same one:
#
#   frame_*  the precursor plus 20 nt on each side ("Precursor_20nt_*").  It is
#            wide enough that a mature shifted by 1 nt still lands inside it,
#            so candidates are built and resolved here.
#   pre_*    the real hairpin (mature + loop + star).  This is what the mature
#            has to fit inside, and the only frame in which the arm-boundary
#            test is meaningful.
#
# Both layouts also carry the mature interval at the same indices.
#   22-col (post-rearrange):     strand 5, frame 6..7 seq 8,  precursor 10..11 seq 12
#   24-col (primary-basic-info): strand 7, frame 5..6 seq 8,  precursor 12..13 seq 10
LAYOUTS = {
    22: {'strand': 5, 'frame_beg': 6, 'frame_end': 7, 'frame_seq': 8,
         'pre_beg': 10, 'pre_end': 11, 'pre_seq': 12},
    24: {'strand': 7, 'frame_beg': 5, 'frame_end': 6, 'frame_seq': 8,
         'pre_beg': 12, 'pre_end': 13, 'pre_seq': 10},
}
MAT_START, MAT_END, MAT_SEQ = 15, 16, 17


def collect_offsets(aln_files, max_3p_gap):
    """
    Return ``{mir_name: {'offsets': Counter, 'seqs': {offset: Counter}}}``.

    ``max_3p_gap`` is the largest 3' shortfall a read may have and still be
    used: a read that stops early says nothing about the 5' boundary.
    """
    offsets = defaultdict(Counter)
    seqs = defaultdict(lambda: defaultdict(Counter))
    for path in aln_files:
        with open(path) as fh:
            for line in fh:
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 5:
                    continue
                match = re.match(r'^(\d+)=(.*)$', fields[2])
                if not match:
                    continue
                ref_len = int(match.group(1))
                name = match.group(2).split('=')[0]
                try:
                    offset = int(fields[3])
                except ValueError:
                    continue
                gap_3p = ref_len - (offset + len(fields[4]))
                if gap_3p > max_3p_gap:
                    continue
                read_weight = re.search(r'_x([\d.]+)', fields[0])
                weight = float(read_weight.group(1)) if read_weight else 1.0
                offsets[name][offset] += weight
                seqs[name][offset][fields[4].upper()] += weight
    return offsets, seqs


def resolve_aln_files(aln_arg, aln_dir):
    files = []
    if aln_arg:
        files.extend(p.strip() for p in aln_arg.split(',') if p.strip())
    if aln_dir:
        files.extend(sorted(glob.glob(os.path.join(aln_dir, '*.mature.aln'))))
    return [f for f in files if os.path.isfile(f)]


def precursor_bounds(parts, beg_key, end_key, seq):
    """
    Read a window's bounds and verify them against the window's sequence.

    basic-info uses a 1-based start with an exclusive end, so the span equals
    the sequence length.  When it does not, the table is not the layout this
    module assumes and the run stops instead of silently shifting a frame.
    """
    try:
        beg = int(parts[beg_key])
        end = int(parts[end_key])
    except (ValueError, IndexError):
        return None
    if beg <= 0 or end <= beg:
        return None
    if end - beg != len(seq):
        return None
    return beg, end


def audit(basic_info, aln_files, *, min_weight, min_fraction, max_shift,
          max_3p_gap, want_corrected, allow_extend=False):
    offsets, seqs = collect_offsets(aln_files, max_3p_gap)
    report = []
    summary = Counter()
    overflow = []
    out_lines = []
    changed = 0

    with open(basic_info) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            layout = LAYOUTS.get(len(parts))
            if layout is None:
                sys.exit(f"Error: {basic_info} has {len(parts)} columns; "
                         f"expected 22 or 24.")
            name = parts[0]
            strand = parts[layout['strand']]
            frame_seq = parts[layout['frame_seq']].upper().replace('U', 'T')
            pre_seq = parts[layout['pre_seq']].upper().replace('U', 'T')
            try:
                mature_start = int(parts[MAT_START])
                mature_end = int(parts[MAT_END])
            except ValueError:
                sys.exit(f"Error: non-integer mature coordinates for {name}")

            frame = precursor_bounds(parts, layout['frame_beg'],
                                     layout['frame_end'], frame_seq)
            if frame is None:
                sys.exit(f"Error: extended-window bounds for {name} do not match "
                         f"its sequence length; refusing to shift a frame.")
            frame_start, frame_end = frame

            pbounds = precursor_bounds(parts, layout['pre_beg'],
                                       layout['pre_end'], pre_seq)
            if pbounds is None:
                sys.exit(f"Error: precursor bounds for {name} do not match its "
                         f"sequence length; refusing to shift a frame.")
            pre_start, pre_end = pbounds

            # Candidates are expressed in the extended window (always wide
            # enough); the arm test and the containment check use the hairpin.
            f_rel_beg, f_rel_end = mature_rel_bounds(frame_start, frame_end,
                                                     mature_start, mature_end,
                                                     strand)
            p_rel_beg, p_rel_end = mature_rel_bounds(pre_start, pre_end,
                                                     mature_start, mature_end,
                                                     strand)
            arm = at_arm_boundary(p_rel_beg, p_rel_end, len(pre_seq))

            verdict = classify(weights=offsets.get(name, Counter()),
                               min_weight=min_weight,
                               min_fraction=min_fraction,
                               max_shift=max_shift,
                               allow_extend=allow_extend)
            status = verdict['status']
            action = ''
            new_start, new_end = mature_start, mature_end
            new_seq = parts[MAT_SEQ]

            if status == SHIFT:
                # The heaviest read at the dominant offset is the witness: the
                # corrected interval has to reproduce it exactly.
                offset = verdict['offset']
                dominant = seqs.get(name, {}).get(offset)
                read_seq = dominant.most_common(1)[0][0] if dominant else ''
                resolved = resolve_candidate(pre_seq=frame_seq,
                                             rel_beg=f_rel_beg,
                                             rel_end=f_rel_end,
                                             shift=verdict['shift'],
                                             read_seq=read_seq)
                coords = None
                if resolved is not None:
                    action, nbeg, nend, nseq = resolved
                    coords = corrected_genomic(pre_start=frame_start,
                                               pre_end=frame_end, strand=strand,
                                               rel_beg=nbeg, rel_end=nend)
                if coords is None:
                    status = OUT_OF_RANGE
                    action = ''
                else:
                    new_start, new_end = coords
                    new_seq = nseq
                    # The corrected mature must still fit inside the annotated
                    # hairpin.  It does not when the precursor was excised flush
                    # with the old (wrong) 3' end, so the precursor end is the
                    # thing that is short.  Report it rather than hide it.
                    if new_start < pre_start or new_end > pre_end:
                        overflow.append(name)
                        summary['precursor_overflow'] += 1

            summary[status] += 1
            report.append([
                name, status,
                '' if verdict['offset'] is None else str(verdict['offset']),
                f"{verdict['fraction']:.3f}",
                f"{verdict['total']:.1f}",
                arm, action,
                str(mature_start), str(mature_end), parts[MAT_SEQ],
                str(new_start), str(new_end), new_seq,
            ])

            if want_corrected and (new_start, new_end, new_seq) != (
                    mature_start, mature_end, parts[MAT_SEQ]):
                parts[MAT_START] = str(new_start)
                parts[MAT_END] = str(new_end)
                parts[MAT_SEQ] = new_seq
                changed += 1
            out_lines.append('\t'.join(parts) + '\n')

    return report, summary, out_lines, changed, overflow


def main():
    parser = argparse.ArgumentParser(
        description="Audit the mature 5' boundary of a basic-info table "
                    "against its bowtie alignment.")
    parser.add_argument('-i', '--input', required=True,
                        help="basic-info table (22 or 24 columns)")
    parser.add_argument('-aln', '--alignment', default=None,
                        help="comma separated *.mature.aln files")
    parser.add_argument('-aln-dir', '--alignment-dir', default=None,
                        help="directory containing *.mature.aln files")
    parser.add_argument('-o', '--output', required=True,
                        help="audit report (TSV)")
    parser.add_argument('-oc', '--output-basic-info', default=None,
                        help="corrected basic-info (requires --apply)")
    parser.add_argument('--apply', action='store_true',
                        help="write the corrected basic-info instead of only "
                             "reporting")
    parser.add_argument('--min-weight', type=float, default=10.0,
                        help="minimum read weight for a locus to be judged "
                             "(default: %(default)s)")
    parser.add_argument('--min-fraction', type=float, default=0.60,
                        help="minimum share of the dominant offset "
                             "(default: %(default)s)")
    parser.add_argument('--max-shift', type=int, default=MAX_SHIFT,
                        help="largest correction in bp (default: %(default)s)")
    parser.add_argument('--allow-extend', action='store_true',
                        help="also extend the mature 5' end when the reads "
                             "start upstream of it (off by default: such a "
                             "locus is reported, not rewritten)")
    parser.add_argument('--max-3p-gap', type=int, default=2,
                        help="largest 3' shortfall of a read still used for the "
                             "5' estimate (default: %(default)s)")
    args = parser.parse_args()

    if args.apply and not args.output_basic_info:
        sys.exit("Error: --apply needs -oc/--output-basic-info.")
    aln_files = resolve_aln_files(args.alignment, args.alignment_dir)
    if not aln_files:
        sys.exit("Error: no *.mature.aln found; pass -aln or -aln-dir.")

    report, summary, out_lines, changed, overflow = audit(
        args.input, aln_files,
        min_weight=args.min_weight,
        min_fraction=args.min_fraction,
        max_shift=args.max_shift,
        max_3p_gap=args.max_3p_gap,
        want_corrected=args.apply,
        allow_extend=args.allow_extend,
    )

    with open(args.output, 'w') as fh:
        fh.write('\t'.join(REPORT_HEADER) + '\n')
        for row in report:
            fh.write('\t'.join(row) + '\n')

    total = len(report)
    flagged = total - summary[OK] - summary[LOW_COVERAGE]
    print(f"[boundary] alignments used: {len(aln_files)}", file=sys.stderr)
    print(f"[boundary] loci: {total}  canonical: {summary[OK]}  "
          f"low coverage: {summary[LOW_COVERAGE]}  flagged: {flagged}",
          file=sys.stderr)
    for status in (SHIFT, SHIFT_UPSTREAM, SHIFT_BEYOND_LIMIT, AMBIGUOUS,
                   OUT_OF_RANGE):
        if summary[status]:
            print(f"[boundary]   {status}: {summary[status]}", file=sys.stderr)
    if overflow:
        print(f"[boundary]   corrected mature overhangs the annotated precursor: "
              f"{len(overflow)} ({', '.join(sorted(overflow))}); the precursor "
              f"end is short by the same amount", file=sys.stderr)

    if args.apply:
        if changed:
            with open(args.output_basic_info, 'w') as fh:
                fh.writelines(out_lines)
            print(f"[boundary] corrected {changed} mature intervals -> "
                  f"{args.output_basic_info}", file=sys.stderr)
        else:
            print("[boundary] no correction needed; basic-info left as is",
                  file=sys.stderr)
        print(f"[boundary] report: {args.output}", file=sys.stderr)
    else:
        print(f"[boundary] report only (use --apply to rewrite). report: "
              f"{args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
