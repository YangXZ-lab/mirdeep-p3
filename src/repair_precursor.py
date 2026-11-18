#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""Re-extend precursors that no longer contain their (corrected) mature.

Reads a 22-column ``*-basic-info`` table, finds every locus whose mature
interval is not fully inside its precursor interval, widens the precursor from
the genome, re-folds the precursor and its +/-20 nt window with RNAfold,
recomputes the star, and writes the repaired table.

Only the columns that describe those two windows and the star are rewritten:
  * the +/-20 nt stem window  (sequence and structure)
  * the precursor             (bounds, sequence and structure)
  * the star                  (sequence and bounds, only if it moved)
The mature interval, the mature name, the family and the MIR_ID are never
touched.

Inputs:
  -i:   basic-info table (22 columns)
  -g:   genome FASTA (random access, no .fai needed)
  -o:   repair report (TSV)
  -oc:  repaired basic-info (written only with --apply)

Example:
  repair_precursor.py -i case1-basic-info -g IRGSP-1.0_genome.fasta \\
      -o case1-precursor-repair.tsv --apply -oc case1-basic-info.repaired
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils.coords import compute_star  # noqa: E402
from utils.genome_fasta import GenomeFasta  # noqa: E402
from utils.precursor_repair import (  # noqa: E402
    CONTAINED, EXTENDED, MAX_EXTENSION, NO_FOLD, SEQ_MISMATCH, STAR_UNRESOLVED,
    TOO_LONG, apply_repair, fold, repair_locus,
)

# 0-based column indices of the 22-column layout produced by annotation.
LAYOUT22 = {
    'chrom': 4, 'strand': 5,
    'stem_beg': 6, 'stem_end': 7, 'stem_seq': 8, 'stem_struc': 9,
    'pre_beg': 10, 'pre_end': 11, 'pre_seq': 12, 'pre_struc': 13,
    'mat_name': 14, 'mat_beg': 15, 'mat_end': 16, 'mat_seq': 17,
    'star_name': 18, 'star_beg': 19, 'star_end': 20, 'star_seq': 21,
}

REPORT_HEADER = [
    'miRNA', 'status', 'old_pre_start', 'old_pre_end', 'old_pre_len',
    'new_pre_start', 'new_pre_end', 'new_pre_len', 'added_low', 'added_high',
    'mature_start', 'mature_end', 'star_start', 'star_end', 'star_moved',
    'columns_changed', 'note',
]


def main():
    parser = argparse.ArgumentParser(
        description="Re-extend precursors that no longer contain their mature "
                    "and re-fold them with RNAfold.")
    parser.add_argument('-i', '--input', required=True,
                        help="basic-info table (22 columns)")
    parser.add_argument('-g', '--genome', required=True,
                        help="genome FASTA used by the annotation run")
    parser.add_argument('-o', '--output', required=True,
                        help="repair report (TSV)")
    parser.add_argument('-oc', '--output-basic-info', default=None,
                        help="repaired basic-info (requires --apply)")
    parser.add_argument('--apply', action='store_true',
                        help="write the repaired table instead of only reporting")
    parser.add_argument('--max-extension', type=int, default=MAX_EXTENSION,
                        help="largest widening allowed on one side "
                             "(default: %(default)s)")
    parser.add_argument('--RNAfold', default='RNAfold',
                        help="RNAfold executable (default: %(default)s)")
    args = parser.parse_args()

    if args.apply and not args.output_basic_info:
        sys.exit("Error: --apply needs -oc/--output-basic-info.")

    def fold_fn(seq):
        return fold(seq, args.RNAfold)

    rows = []
    with open(args.input) as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 22:
                sys.exit(f"Error: {args.input} line {lineno} has {len(parts)} "
                         f"columns; this tool repairs the 22-column layout only.")
            rows.append(parts)

    report = []
    out_lines = []
    counts = {}
    total_changed = 0

    with GenomeFasta(args.genome) as genome:
        for parts in rows:
            name = parts[0]
            chrom = parts[LAYOUT22['chrom']]
            strand = parts[LAYOUT22['strand']]
            pre_start = int(parts[LAYOUT22['pre_beg']])
            pre_end = int(parts[LAYOUT22['pre_end']])
            mat_start = int(parts[LAYOUT22['mat_beg']])
            mat_end = int(parts[LAYOUT22['mat_end']])
            old_pre_seq = parts[LAYOUT22['pre_seq']].upper()
            mature_seq = parts[LAYOUT22['mat_seq']].upper()

            if chrom not in genome:
                counts['chrom_missing'] = counts.get('chrom_missing', 0) + 1
                report.append([name, 'chrom_missing', pre_start, pre_end,
                               len(old_pre_seq), '', '', '', 0, 0,
                               mat_start, mat_end, '', '', '', 0,
                               f"{chrom} not in the genome FASTA"])
                out_lines.append('\t'.join(parts) + '\n')
                continue

            chrom_len = genome.length(chrom)

            def fetch(start, end, st, _c=chrom):
                return genome.fetch(_c, start, end, st)

            plan = repair_locus(
                pre_start=pre_start, pre_end=pre_end,
                mature_start=mat_start, mature_end=mat_end,
                mature_seq=mature_seq, strand=strand,
                old_pre_seq=old_pre_seq,
                old_stem_seq=parts[LAYOUT22['stem_seq']],
                old_stem_struc=parts[LAYOUT22['stem_struc']],
                fetch=fetch, chrom_len=chrom_len,
                compute_star=compute_star, fold_fn=fold_fn,
                max_extension=args.max_extension,
            )

            status = plan['status']
            counts[status] = counts.get(status, 0) + 1
            changed = 0
            if args.apply and status == EXTENDED:
                changed = apply_repair(parts, plan, LAYOUT22)
                total_changed += 1

            old_star_beg = parts[LAYOUT22['star_beg']]
            old_star_end = parts[LAYOUT22['star_end']]
            new_star_beg = plan['new_star_start']
            new_star_end = plan['new_star_end']
            if new_star_beg is None:
                star_moved = ''
            else:
                star_moved = ('no' if (str(new_star_beg) == old_star_beg
                                       and str(new_star_end) == old_star_end)
                              else 'yes')
            report.append([
                name, status,
                pre_start, pre_end, len(old_pre_seq),
                plan['new_pre_start'], plan['new_pre_end'],
                len(plan['new_pre_seq'] or old_pre_seq),
                plan['added_low'], plan['added_high'],
                mat_start, mat_end,
                new_star_beg if new_star_beg is not None else '',
                new_star_end if new_star_end is not None else '',
                star_moved, changed, plan['note'],
            ])
            out_lines.append('\t'.join(parts) + '\n')

    with open(args.output, 'w') as fh:
        fh.write('\t'.join(REPORT_HEADER) + '\n')
        for row in report:
            fh.write('\t'.join(str(x) for x in row) + '\n')

    total = len(report)
    overhang = total - counts.get(CONTAINED, 0) - counts.get('chrom_missing', 0)
    print(f"[precursor] loci: {total}  mature inside precursor: "
          f"{counts.get(CONTAINED, 0)}  needing repair: {overhang}",
          file=sys.stderr)
    for status in (EXTENDED, TOO_LONG, SEQ_MISMATCH, STAR_UNRESOLVED, NO_FOLD):
        if counts.get(status):
            print(f"[precursor]   {status}: {counts[status]}", file=sys.stderr)
    if counts.get('chrom_missing'):
        print(f"[precursor]   chrom_missing: {counts['chrom_missing']}",
              file=sys.stderr)
    if args.apply:
        if total_changed:
            with open(args.output_basic_info, 'w') as fh:
                fh.writelines(out_lines)
            print(f"[precursor] repaired {total_changed} precursor/star "
                  f"records -> {args.output_basic_info}", file=sys.stderr)
        else:
            print("[precursor] nothing to repair; table left as is",
                  file=sys.stderr)
    print(f"[precursor] report: {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
