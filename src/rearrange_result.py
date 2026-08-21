#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
Rearrange columns of a miRNA primary information file.

"""

import sys
import argparse
from typing import List, Optional


def reorder_line(fields: List[str]) -> Optional[str]:
    """
    Reorder fields according to the original Perl script.

    Input fields indices (0-based, expected 24 columns):
        0: miRNA ID (e.g., N710-MIR156a)
        1: MIR_ID (e.g., MIR_ID000000001)
        2: MIR family (e.g., MIR156)
        3: Species
        4: Chromosome
        5: Pre_start_1 (original position)
        6: Pre_end_1
        7: Strand
        8: Long precursor sequence (with flanking)
        9: Long secondary structure
        10: Short sequence (mature-containing region)
        11: Short secondary structure
        12: Coordinate_1
        13: Coordinate_2
        14: Mature ID
        15: Mature start
        16: Mature end
        17: Mature sequence
        18: Star ID
        19: Star sequence
        20: Star start
        21: Star end
        22: Confidence (or other, skipped)
        23: Source (skipped)

    Output order (22 columns):
        0,1,2,3,4,7,5,6,8,9,12,13,10,11,14,15,16,17,18,20,21,19

    Args:
        fields: List of fields from a tab-separated line.

    Returns:
        Tab-separated reordered line, or None if insufficient fields.
    """
    if len(fields) < 24:
        print(f"Warning: Line has {len(fields)} fields, expected 24. Skipping.",
              file=sys.stderr)
        return None

    # Define output order indices (0-based)
    order = [0, 1, 2, 3, 4, 7, 5, 6, 8, 9, 12, 13, 10, 11, 14, 15, 16, 17, 18, 20, 21, 19]

    out_fields = [fields[i] for i in order]
    return '\t'.join(out_fields)


def process_file(input_file: str, output_file: str) -> None:
    """
    Read input file, reorder each line, and write to output file.

    Args:
        input_file: Path to input file.
        output_file: Path to output file.
    """
    processed = 0
    skipped = 0

    try:
        with open(input_file, 'r') as infh, open(output_file, 'w') as outfh:
            for line_num, line in enumerate(infh, 1):
                line = line.rstrip('\n')
                if not line:
                    continue

                fields = line.split('\t')
                out_line = reorder_line(fields)
                if out_line is None:
                    skipped += 1
                    continue

                outfh.write(out_line + '\n')
                processed += 1

    except Exception as e:
        print(f"Error: Failed to process files: {e}", file=sys.stderr)
        sys.exit(1)

    # Summary
    print(f"Processing completed.", file=sys.stderr)
    print(f"  Lines processed: {processed}", file=sys.stderr)
    if skipped:
        print(f"  Lines skipped: {skipped}", file=sys.stderr)
    print(f"  Output written to: {output_file}", file=sys.stderr)


def main():
    """Command line interface."""
    parser = argparse.ArgumentParser(
        description='Rearrange columns of a miRNA primary information file.',
        epilog='Input: 24-column tab-separated file.\n'
               'Output: 22-column tab-separated file with columns reordered.',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input file path (24 columns).'
    )

    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output file path (reordered 22 columns).'
    )

    args = parser.parse_args()

    process_file(args.input, args.output)


if __name__ == "__main__":
    main()
