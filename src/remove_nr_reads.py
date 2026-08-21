#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
Select non-redundant predictions based on the highest expression value.

"""

import sys
import re
import argparse
from typing import Dict, List, Set


def extract_numeric_value(id_str: str) -> float:
    """
    Extract the numeric value from an ID string.

    Expected pattern: '_x' followed by a number (integer or decimal),
    then a hyphen '-'. Example: "N710-Chr01_76-read01601981_x1443.6683959533646-4"
    returns 1443.6683959533646.

    Args:
        id_str: The ID string.

    Returns:
        The extracted float value. If no match is found, returns -inf
        so that such IDs are never selected as the highest.

    Raises:
        ValueError: If the pattern is not found (handled by returning -inf).
    """
    match = re.search(r'_x(\d+(?:\.\d+)?)-', id_str)
    if match:
        return float(match.group(1))
    # If no pattern found, return negative infinity so it won't be chosen
    return float('-inf')


def parse_merge_file(merge_file: str) -> Dict[str, str]:
    """
    Parse the merge_output.bed file to determine the representative ID for each group.

    Args:
        merge_file: Path to the merge_output.bed file.

    Returns:
        A dictionary mapping each selected representative ID to a marker ("T").
        The marker is not used further, but kept for compatibility with the Perl version.
    """
    representatives: Dict[str, str] = {}
    line_count = 0
    skipped_lines = 0

    try:
        with open(merge_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.rstrip('\n')
                if not line:
                    continue

                fields = line.split('\t')
                if len(fields) < 4:
                    print(f"Warning: Skipping line {line_num} in {merge_file}: "
                          f"expected at least 4 columns, got {len(fields)}",
                          file=sys.stderr)
                    skipped_lines += 1
                    continue

                id_list_str = fields[3]  # 4th column (0-index)
                id_list = [id_str.strip() for id_str in id_list_str.split(',') if id_str.strip()]

                if not id_list:
                    print(f"Warning: Empty ID list at line {line_num}, skipping.",
                          file=sys.stderr)
                    skipped_lines += 1
                    continue

                highest_value = float('-inf')
                best_id = None

                for id_str in id_list:
                    value = extract_numeric_value(id_str)
                    if value > highest_value:
                        highest_value = value
                        best_id = id_str

                if best_id is not None:
                    representatives[best_id] = "T"
                    line_count += 1
                else:
                    print(f"Warning: Could not extract any valid numeric value "
                          f"from ID list at line {line_num}, skipping.",
                          file=sys.stderr)
                    skipped_lines += 1

    except Exception as e:
        print(f"Error: Cannot read merge file {merge_file}: {e}", file=sys.stderr)
        sys.exit(1)

    if line_count == 0:
        print(f"Error: No valid groups found in {merge_file}", file=sys.stderr)
        sys.exit(1)

    if skipped_lines:
        print(f"  (Skipped {skipped_lines} malformed lines in merge file)", file=sys.stderr)

    return representatives


def filter_predictions(pred_file: str, representatives: Dict[str, str], output_file: str) -> None:
    """
    Filter the all-mod_fp_prediction file, keeping only rows whose ID is a representative.

    Args:
        pred_file: Path to the prediction file (8 columns, 4th column is ID).
        representatives: Dictionary of representative IDs (keys) with any value.
        output_file: Path to the output file.
    """
    kept_count = 0
    total_count = 0
    rep_set = set(representatives.keys())

    try:
        with open(pred_file, 'r') as infh, open(output_file, 'w') as outfh:
            for line in infh:
                line = line.rstrip('\n')
                if not line:
                    continue
                total_count += 1

                fields = line.split('\t')
                if len(fields) < 4:
                    print(f"Warning: Skipping line with <4 columns in {pred_file}",
                          file=sys.stderr)
                    continue

                id_field = fields[3]  # 4th column
                if id_field in rep_set:
                    outfh.write(line + '\n')
                    kept_count += 1

    except Exception as e:
        print(f"Error: Cannot process prediction file {pred_file}: {e}", file=sys.stderr)
        sys.exit(1)

    # Summary
    print(f"\nFiltering completed.", file=sys.stderr)
    print(f"  Total lines in prediction file: {total_count}", file=sys.stderr)
    print(f"  Lines kept (representatives): {kept_count}", file=sys.stderr)
    print(f"  Output written to: {output_file}", file=sys.stderr)


def main():
    """Command line interface."""
    parser = argparse.ArgumentParser(
        description='Select non-redundant predictions based on highest expression value.',
        epilog='For each group in the merge file (comma-separated IDs in column 4), '
               'the ID with the highest value (extracted from "_x<number>-") is chosen. '
               'Then only those rows from the prediction file are output.',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '-m', '--merge',
        required=True,
        help='Merge output file (BED-like, 6 columns, 4th column is comma-separated IDs).'
    )

    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Prediction file (all-mod_fp_prediction, 8 columns, 4th column is ID).'
    )

    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output file (non-redundant predictions, same format as input).'
    )

    args = parser.parse_args()

    # Parse merge file to get representatives
    print(f"Reading merge file: {args.merge}", file=sys.stderr)
    representatives = parse_merge_file(args.merge)
    print(f"  Found {len(representatives)} representative IDs.", file=sys.stderr)

    # Filter prediction file
    print(f"\nFiltering prediction file: {args.input}", file=sys.stderr)
    filter_predictions(args.input, representatives, args.output)


if __name__ == "__main__":
    main()
