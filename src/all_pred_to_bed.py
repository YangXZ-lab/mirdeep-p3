#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
Convert P_prediction (or all-mod_fp_prediction) files to a BED-like format.

"""

import sys
import re
import argparse
from typing import Tuple, Optional


def extract_start_from_range(range_str: str) -> Optional[int]:
    """
    Extract the start coordinate from a range string formatted as "start..end".

    Args:
        range_str: String like "222475172..222475192"

    Returns:
        Start coordinate as integer, or None if format is invalid.
    """
    match = re.match(r'^(\d+)\.\.(\d+)$', range_str)
    if match:
        return int(match.group(1))
    return None


def determine_signature(mature_beg: int, mature_end: int,
                        pre_beg: int, pre_end: int,
                        strand: str) -> str:
    """
    Determine the position signature of mature miRNA relative to pre-miRNA.

    Args:
        mature_beg: Adjusted start coordinate of mature miRNA.
        mature_end: Adjusted end coordinate of mature miRNA.
        pre_beg: Adjusted start coordinate of pre-miRNA.
        pre_end: Adjusted end coordinate of pre-miRNA.
        strand: Strand direction ('+' or '-').

    Returns:
        '5' if mature aligns with 5' end of pre-miRNA,
        '3' if aligns with 3' end,
        'A' for ambiguous (internal) location.
    """
    if (mature_beg == pre_beg and strand == '+') or \
       (mature_end == pre_end and strand == '-'):
        return '5'
    elif (mature_beg == pre_beg and strand == '-') or \
         (mature_end == pre_end and strand == '+'):
        return '3'
    else:
        return 'A'


def process_line(fields: list) -> Optional[str]:
    """
    Process a single line from the input prediction file.

    Args:
        fields: List of 8 tab-separated fields.

    Returns:
        Formatted output line (6 columns, tab-separated) or None if line is invalid.
    """
    if len(fields) < 8:
        print(f"Warning: Skipping line with {len(fields)} columns (expected 8): {fields}",
              file=sys.stderr)
        return None

    chrom = fields[0]
    strand = fields[1]
    name = fields[3]          # identifier
    mature_range = fields[4]  # e.g., "222475172..222475192"
    pre_range = fields[5]     # e.g., "222475117..222475192"
    mature_seq = fields[6]
    pre_seq = fields[7]

    # Extract start coordinates from ranges
    mature_beg = extract_start_from_range(mature_range)
    pre_beg = extract_start_from_range(pre_range)

    if mature_beg is None or pre_beg is None:
        print(f"Warning: Cannot parse coordinate range(s): mature='{mature_range}', pre='{pre_range}'",
              file=sys.stderr)
        return None

    # Adjust coordinates based on strand
    if strand == '+':
        pre_beg -= 2
        mature_beg -= 2
    elif strand == '-':
        pre_beg += 1
        mature_beg += 1
    else:
        print(f"Warning: Unknown strand '{strand}', skipping line.", file=sys.stderr)
        return None

    # Calculate lengths and end coordinates
    mature_len = len(mature_seq)
    pre_len = len(pre_seq)
    mature_end = mature_beg + mature_len
    pre_end = pre_beg + pre_len

    # Determine signature (5'/3'/A)
    sign = determine_signature(mature_beg, mature_end, pre_beg, pre_end, strand)

    # Output: chrom, pre_start, pre_end, name, sign, strand
    out_line = f"{chrom}\t{pre_beg}\t{pre_end}\t{name}\t{sign}\t{strand}"
    return out_line


def convert_predictions_to_bed(input_file: str, output_file: str) -> None:
    """
    Convert the prediction file to BED-like format.

    Args:
        input_file: Path to input file (8-column tab-separated).
        output_file: Path to output file (6-column BED-like).
    """
    processed_lines = 0
    skipped_lines = 0

    try:
        with open(input_file, 'r') as infh, open(output_file, 'w') as outfh:
            for line_num, line in enumerate(infh, 1):
                line = line.rstrip('\n')
                if not line:
                    continue

                fields = line.split('\t')
                out_line = process_line(fields)

                if out_line is None:
                    skipped_lines += 1
                    continue

                outfh.write(out_line + '\n')
                processed_lines += 1

    except Exception as e:
        print(f"Error: Failed to process file {input_file}: {e}", file=sys.stderr)
        sys.exit(1)

    # Print summary to stderr
    print(f"Conversion completed.", file=sys.stderr)
    print(f"  Processed lines: {processed_lines}", file=sys.stderr)
    if skipped_lines:
        print(f"  Skipped lines: {skipped_lines}", file=sys.stderr)
    print(f"  Output written to: {output_file}", file=sys.stderr)


def main():
    """Command line interface."""
    parser = argparse.ArgumentParser(
        description='Convert P_prediction (or all-mod_fp_prediction) files to BED-like format.',
        epilog='Output is 6-column: chrom<TAB>pre_start<TAB>pre_end<TAB>name<TAB>sign<TAB>strand',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input file path (8-column tab-separated prediction file).'
    )

    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output file path (BED-like format).'
    )

    args = parser.parse_args()

    convert_predictions_to_bed(args.input, args.output)


if __name__ == "__main__":
    main()
