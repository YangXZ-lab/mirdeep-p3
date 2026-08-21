#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
Sort and deduplicate a BED-like file by the 4th column (identifier).

"""

import sys
import argparse
from typing import Dict, List, Tuple


def parse_bed_line(line: str) -> Tuple[str, int, int, str, str]:
    """
    Parse a line of the input BED-like file.

    Expected format: chrom<TAB>start<TAB>end<TAB>name<TAB>sign<TAB>strand

    Args:
        line: A line from the input file.

    Returns:
        A tuple (chrom, start, end, name, strand).

    Raises:
        ValueError: If the line does not have exactly 6 fields or
                    start/end cannot be converted to integers.
    """
    fields = line.strip().split('\t')
    if len(fields) != 6:
        raise ValueError(f"Expected 6 columns, got {len(fields)}")

    chrom = fields[0]
    try:
        start = int(fields[1])
        end = int(fields[2])
    except ValueError as e:
        raise ValueError(f"Invalid start/end values: {fields[1]}, {fields[2]}") from e

    name = fields[3]
    # sign = fields[4]  # not used in output
    strand = fields[5]

    return chrom, start, end, name, strand


def process_bed_file(input_file: str, output_file: str) -> None:
    """
    Read, deduplicate, sort, and write the BED-like file.

    Args:
        input_file: Path to the input file (6 columns).
        output_file: Path to the output file.
    """
    records: Dict[str, Tuple[str, int, int, str]] = {}
    line_count = 0
    dup_count = 0

    try:
        with open(input_file, 'r') as infh:
            for line_num, line in enumerate(infh, 1):
                line = line.rstrip('\n')
                if not line:
                    continue

                try:
                    chrom, start, end, name, strand = parse_bed_line(line)
                except ValueError as e:
                    print(f"Warning: Skipping line {line_num}: {e}", file=sys.stderr)
                    continue

                if name in records:
                    dup_count += 1
                records[name] = (chrom, start, end, strand)
                line_count += 1

    except Exception as e:
        print(f"Error: Cannot read input file {input_file}: {e}", file=sys.stderr)
        sys.exit(1)

    if not records:
        print("Error: No valid records found in input file.", file=sys.stderr)
        sys.exit(1)

    # Sort by: chromosome (lexicographically), start (numerical), end (numerical)
    sorted_names = sorted(
        records.keys(),
        key=lambda name: (records[name][0], records[name][1], records[name][2])
    )

    try:
        with open(output_file, 'w') as outfh:
            for name in sorted_names:
                chrom, start, end, strand = records[name]
                # Output: chrom, start, end, name, '.', strand
                outfh.write(f"{chrom}\t{start}\t{end}\t{name}\t.\t{strand}\n")
    except Exception as e:
        print(f"Error: Cannot write output file {output_file}: {e}", file=sys.stderr)
        sys.exit(1)

    # Summary
    print(f"Processing completed.", file=sys.stderr)
    print(f"  Total lines read (valid): {line_count}", file=sys.stderr)
    if dup_count:
        print(f"  Duplicate names removed: {dup_count}", file=sys.stderr)
    print(f"  Unique records written: {len(records)}", file=sys.stderr)
    print(f"  Output written to: {output_file}", file=sys.stderr)


def main():
    """Command line interface."""
    parser = argparse.ArgumentParser(
        description='Sort and deduplicate a BED-like file by the 4th column (name).',
        epilog='Input: 6-column tab-separated (chrom, start, end, name, sign, strand)\n'
               'Output: same format but with 5th column replaced by "." and sorted.',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input BED-like file path (6 columns).'
    )

    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output BED-like file path (sorted and deduplicated).'
    )

    args = parser.parse_args()

    process_bed_file(args.input, args.output)


if __name__ == "__main__":
    main()
