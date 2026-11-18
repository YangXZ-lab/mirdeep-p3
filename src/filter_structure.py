#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""Filter out loci whose mature and star cannot both sit on the precursor.

Reads a 22-column ``*-basic-info`` table and drops every row whose two arms
overlap, whose loop between them is shorter than ``--min-arm-gap``, or whose
structure is missing.  Writes the filtered table and a report naming each
rejected locus and why, so a shrinking table is always explainable.

This runs right after basic-info is produced, before any alignment,
quantification or clustering, so no work is spent on a locus that will not
survive.

Inputs:
  -i:   basic-info table (22 columns)
  -o:   filter report (TSV), one row per input locus
  -oc:  filtered basic-info (written only with --apply)

Example:
  filter_structure.py -i case1-basic-info -o case1-structure-filter.tsv \\
      --apply -oc case1-basic-info.filtered
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils.structure_filter import (  # noqa: E402
    ARM_OVERLAP, BAD_COORDINATES, KEPT, MIN_ARM_GAP, MISSING_PRECURSOR_STRUCTURE,
    MISSING_STAR, MISSING_STEM_STRUCTURE, REPORT_HEADER, SHORT_ARM_GAP,
    filter_rows, report_line, summarize,
)

# The report lists the interesting statuses first, then the noise.
_STATUS_ORDER = (ARM_OVERLAP, SHORT_ARM_GAP, MISSING_PRECURSOR_STRUCTURE,
                 MISSING_STEM_STRUCTURE, MISSING_STAR, BAD_COORDINATES)


def main():
    parser = argparse.ArgumentParser(
        description="Drop loci whose mature and star overlap, or whose loop "
                    "between the arms is too short.")
    parser.add_argument('-i', '--input', required=True,
                        help="basic-info table (22 columns)")
    parser.add_argument('-o', '--output', required=True,
                        help="filter report (TSV), one row per input locus")
    parser.add_argument('-oc', '--output-basic-info', default=None,
                        help="filtered basic-info (requires --apply)")
    parser.add_argument('--apply', action='store_true',
                        help="write the filtered table instead of only reporting")
    parser.add_argument('--min-arm-gap', type=int, default=MIN_ARM_GAP,
                        help="minimum loop length between mature and star, in "
                             "bases (default: %(default)s). Use 0 to keep every "
                             "non-overlapping locus.")
    args = parser.parse_args()

    if args.min_arm_gap < 0:
        sys.exit("Error: --min-arm-gap must be non-negative.")
    if args.apply and not args.output_basic_info:
        sys.exit("Error: --apply needs -oc/--output-basic-info.")

    rows = []
    with open(args.input) as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 22:
                sys.exit(f"Error: {args.input} line {lineno} has {len(parts)} "
                         f"columns; this filter works on the 22-column layout.")
            rows.append(parts)

    kept, report = filter_rows(rows, min_gap=args.min_arm_gap)
    counts = summarize(report)
    removed = len(rows) - len(kept)

    with open(args.output, 'w') as fh:
        fh.write('\t'.join(REPORT_HEADER) + '\n')
        for entry in report:
            fh.write(report_line(entry) + '\n')

    print(f"[structure] loci: {len(rows)}  kept: {len(kept)}  removed: {removed} "
          f"(min arm gap: {args.min_arm_gap} nt)", file=sys.stderr)
    for status in _STATUS_ORDER:
        if counts.get(status):
            print(f"[structure]   {status}: {counts[status]}", file=sys.stderr)
    if removed:
        dropped = [str(e['miRNA']) for e in report if e['status'] != KEPT]
        shown = ', '.join(dropped[:12])
        if len(dropped) > 12:
            shown += f", ... (+{len(dropped) - 12})"
        print(f"[structure]   removed: {shown}", file=sys.stderr)

    if args.apply:
        if not kept:
            sys.exit("Error: every locus was filtered out; refusing to write an "
                     "empty basic-info.")
        with open(args.output_basic_info, 'w') as fh:
            for parts in kept:
                fh.write('\t'.join(parts) + '\n')
        print(f"[structure] filtered table -> {args.output_basic_info}",
              file=sys.stderr)
    print(f"[structure] report: {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
