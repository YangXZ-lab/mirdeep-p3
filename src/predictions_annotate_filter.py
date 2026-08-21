#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
predictions_annotate_filter.py

Filter and annotate miRDP2 miRNA prediction results (TSV format with 9 columns).
Based on the original Perl one-liner logic but adjusted for actual column count:
  - Column 9 (index 8) contains the conservation flag ("non_conserved" or "conserved=...").
  - Column 7 (index 6): mature miRNA sequence (used for length check 21-22).
  - Column 8 (index 7): precursor sequence (used for length check >=56).
  - Column 3 (index 2): read identifier with "_x<number>" suffix; number must be >=10
    for a specific fallback branch.

Processing rules:
  If column 9 contains "non_conserved":
      If len(col7) in [21, 22]:
          Keep row if len(col8) >= 56 (output unchanged).
      Else (len(col7) not in 21-22):
          Extract the numeric value after "_x" in col3.
          Keep row if value >= 10 and len(col8) >= 56.
          Output original row + "\t3\tmiRDP2".
  Else (col9 does NOT contain "non_conserved"):
      Keep row unconditionally.
      Output original row + "\t3\tmiRDP2=<col4>".

Rows not meeting any keep condition are discarded.

Usage:
    python predictions_annotate_filter.py -i input.tsv -o output.tsv
"""

import sys
import re
import argparse

def extract_x_value(text):
    """
    Extract the numeric value following '_x' in the given string.
    Example: "read01827486_x11.2422510139651" -> 11.2422510139651 (as float)
    Returns None if pattern not found.
    """
    match = re.search(r'_x(\d+\.?\d*)', text)
    if match:
        return float(match.group(1))
    return None

def process_line(line):
    """
    Process a single TSV line (string). Returns a string to be written to output,
    or None if the line should be discarded.
    """
    line = line.rstrip('\n')
    if not line:
        return None

    fields = line.split('\t')
    if len(fields) < 9:
        sys.stderr.write(f"Warning: line has fewer than 9 columns, skipping: {line}\n")
        return None

    col3 = fields[2]   # read identifier with _x value
    col4 = fields[3]   # precursor ID
    col7 = fields[6]   # mature miRNA sequence
    col8 = fields[7]   # precursor sequence
    col9 = fields[8]   # conservation flag

    is_non_conserved = 'non_conserved' in col9

    if is_non_conserved:
        len_col7 = len(col7)
        len_col8 = len(col8)

        # Branch 1: mature length 20-22
        if 20 <= len_col7 <= 22:
            if len_col8 >= 56:
                return line   # output unchanged
            else:
                return None
        else:
            # Branch 2: mature length outside 20-22
            value = extract_x_value(col3)
            if value is not None and value >= 10.0 and len_col8 >= 56:
                return f"{line}\t3\tmiRDP2"
            else:
                return None
    else:
        # Not non_conserved: keep always, append tag
        return f"{line}\t3\tmiRDP2={col4}"

def main():
    parser = argparse.ArgumentParser(
        description="Filter miRDP2 TSV predictions and add annotation columns."
    )
    parser.add_argument('-i', '--input', required=True,
                        help='Input TSV file (9 columns)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output filtered TSV file')
    args = parser.parse_args()

    try:
        with open(args.input, 'r', encoding='utf-8') as infile, \
             open(args.output, 'w', encoding='utf-8') as outfile:

            kept_count = 0
            for line_num, line in enumerate(infile, 1):
                result = process_line(line)
                if result is not None:
                    outfile.write(result + '\n')
                    kept_count += 1

            sys.stderr.write(f"Processing complete. Kept {kept_count} lines.\n")

    except IOError as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
