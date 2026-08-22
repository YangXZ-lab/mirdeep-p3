#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
Extract columns from a tab-delimited file according to case1 and case2 group definitions.
Both groups must be specified, have equal size, and must not overlap.
The output contains column 1 followed by all columns of case1 then all columns of case2.
"""

import argparse
import sys


def parse_columns(col_str: str):
    """Parse comma-separated column numbers into a list of 1-based integers."""
    parts = col_str.split(',')
    cols = []
    for p in parts:
        p = p.strip()
        if not p:
            raise ValueError("Empty column number in specification")
        try:
            col = int(p)
        except ValueError:
            raise ValueError(f"Invalid column number: '{p}'")
        if col < 1:
            raise ValueError(f"Column number must be >= 1: {col}")
        cols.append(col)
    return cols


def validate_groups(case1_str: str, case2_str: str):
    """Validate and return column lists for case1 and case2."""
    if case1_str is None or case2_str is None:
        raise ValueError("Both --case1 and --case2 must be provided")
    cols1 = parse_columns(case1_str)
    cols2 = parse_columns(case2_str)

    if len(cols1) != len(cols2):
        raise ValueError(
            f"Group sizes differ: case1 has {len(cols1)} columns, case2 has {len(cols2)} columns"
        )
    if set(cols1) & set(cols2):
        raise ValueError("case1 and case2 columns must not overlap")
    return cols1, cols2


def extract_columns(input_path: str, output_path: str,
                    cols1: list, cols2: list):
    """Read input file and write selected columns to output file."""
    # Column order: 1, then case1 columns, then case2 columns
    selected = [1] + cols1 + cols2
    max_col = max(selected)

    with open(input_path, 'r') as fin, open(output_path, 'w') as fout:
        for lineno, line in enumerate(fin, 1):
            line = line.rstrip('\n')
            # Preserve empty lines
            if not line:
                fout.write('\n')
                continue

            fields = line.split('\t')
            if len(fields) < max_col:
                raise ValueError(
                    f"Line {lineno} has {len(fields)} columns, "
                    f"but at least {max_col} are required"
                )
            out_fields = [fields[i - 1] for i in selected]
            fout.write('\t'.join(out_fields) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Extract columns from a tabular file by case groups.'
    )
    parser.add_argument('-i', '--input', required=True,
                        help='Input file path (tab-delimited)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output file path')
    parser.add_argument('--case1',
                        help='Comma-separated 1-based column numbers for group 1')
    parser.add_argument('--case2',
                        help='Comma-separated 1-based column numbers for group 2')
    args = parser.parse_args()

    try:
        case1_cols, case2_cols = validate_groups(args.case1, args.case2)
    except ValueError as e:
        sys.stderr.write(f"Argument error: {e}\n")
        sys.exit(1)

    try:
        extract_columns(args.input, args.output, case1_cols, case2_cols)
    except Exception as e:
        sys.stderr.write(f"Processing error: {e}\n")
        sys.exit(1)


if __name__ == '__main__':
    main()
