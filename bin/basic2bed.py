#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
Extract specific columns from a tab-separated file.
Input columns (1-based): 5, 7, 8, 1, 4, 6
Output order: col5, col7, col8, col1, col4, col6

Usage:
    python extract_columns.py -i input.tsv -o output.tsv
"""

import argparse

def main():
    parser = argparse.ArgumentParser(description='Extract columns 5,7,8,1,4,6 from a TSV file.')
    parser.add_argument('-i', '--input', required=True, help='Input TSV file')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    args = parser.parse_args()

    with open(args.input, 'r') as fin, open(args.output, 'w') as fout:
        for line in fin:
            line = line.rstrip('\n')
            if not line:
                continue
            fields = line.split('\t')
            # Check if there are at least 8 columns (1-based index 8)
            if len(fields) < 8:
                # If line has fewer columns, skip or handle? AWK would print empty fields for missing.
                # To mimic AWK, we could pad with empty strings up to 8.
                while len(fields) < 8:
                    fields.append('')
            # Extract columns (0-based indices: 4,6,7,0,3,5)
            selected = [fields[4], fields[6], fields[7], fields[0], fields[3], fields[5]]
            fout.write('\t'.join(selected) + '\n')

if __name__ == '__main__':
    main()
