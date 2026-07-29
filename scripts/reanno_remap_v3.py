#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
Combine cluster, singleton and alignment results into a three-column
miRNA-to-family mapping file.  Third column (ratio) is preserved for
entries from the --multi file; for other sources it is left empty.

Usage:
    python merge_family.py --alncl alncl.txt --single single.txt --multi multi.txt \
                           -s 2001 -o output.txt
"""

import argparse
import sys
from collections import defaultdict

def parse_alncl(filepath):
    """
    Read the alncl file (two columns: mirna, family)
    Return a list of (mirna, family, extra) with extra=''.
    """
    entries = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                sys.stderr.write(f"Warning: skipping malformed alncl line: {line}\n")
                continue
            mirna = parts[0].strip()
            family = parts[1].strip()
            entries.append((mirna, family, ''))
    return entries

def parse_multi(filepath, start_num):
    """
    Read the multi (cluster) file with optional third column.
    Returns a list of (mirna, family, extra) tuples and the next available
    family number.
    """
    clusters = defaultdict(list)        # cluster_id -> list of (mirna, extra)
    cluster_order = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                sys.stderr.write(f"Warning: skipping malformed multi line: {line}\n")
                continue
            mirna = parts[0].strip()
            cluster_id = parts[1].strip()
            extra = parts[2].strip() if len(parts) >= 3 else ''
            if cluster_id not in clusters:
                cluster_order.append(cluster_id)
            clusters[cluster_id].append((mirna, extra))

    # Sort cluster IDs by their numeric part
    def cluster_sort_key(cid):
        import re
        m = re.search(r'\d+', cid)
        return int(m.group()) if m else 0
    cluster_order.sort(key=cluster_sort_key)

    entries = []
    next_num = start_num
    for cid in cluster_order:
        family = f"MIR{next_num}"
        for mirna, extra in clusters[cid]:
            entries.append((mirna, family, extra))
        next_num += 1

    return entries, next_num

def parse_single(filepath, start_num):
    """
    Read the single file and return a list of (mirna, family, '') tuples.
    Each miRNA gets its own sequential family number.
    """
    entries = []
    with open(filepath, 'r') as f:
        current_num = start_num
        for line in f:
            mirna = line.strip()
            if not mirna:
                continue
            family = f"MIR{current_num}"
            entries.append((mirna, family, ''))
            current_num += 1
    return entries

def main():
    parser = argparse.ArgumentParser(
        description='Create miRNA-to-family mapping from alncl, multi, and single files.')
    parser.add_argument('--alncl', default=None,
                        help='Alncl file (two columns: mirna, family)')
    parser.add_argument('--single', default=None,
                        help='Single file (one miRNA per line)')
    parser.add_argument('--multi', default=None,
                        help='Multi (cluster) file (two or three columns: mirna, cluster_id, ratio)')
    parser.add_argument('-s', '--start', type=int, required=True,
                        help='Starting family number for clusters')
    parser.add_argument('-o', '--output', required=True,
                        help='Output mapping file (three columns: mirna, family, extra)')
    args = parser.parse_args()

    all_entries = []

    # 1. Process alncl file
    if args.alncl:
        all_entries.extend(parse_alncl(args.alncl))

    # 2. Process multi file
    next_single_start = args.start
    if args.multi:
        multi_entries, next_single_start = parse_multi(args.multi, args.start)
        all_entries.extend(multi_entries)

    # 3. Process single file
    if args.single:
        single_entries = parse_single(args.single, next_single_start)
        all_entries.extend(single_entries)

    # Write output: three tab-separated columns; extra may be empty
    with open(args.output, 'w') as out:
        for mirna, family, extra in all_entries:
            out.write(f"{mirna}\t{family}\t{extra}\n")

    print(f"Output written to {args.output}", file=sys.stderr)

if __name__ == '__main__':
    main()
