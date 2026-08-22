#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

import argparse
import sys

def read_fai(fai_file):
    chrom_lengths = {}
    with open(fai_file, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                chrom = parts[0]
                length = int(parts[1])
                chrom_lengths[chrom] = length
    return chrom_lengths

def get_upstream_region(chrom, start, end, strand, chrom_lengths, upstream=2000):
    if chrom not in chrom_lengths:
        print(f"Warning: chromosome '{chrom}' not found in .fai file, skipping.", file=sys.stderr)
        return None

    chr_len = chrom_lengths[chrom]

    if strand == '+':
        new_start = max(0, start - upstream)
        new_end = start
    elif strand == '-':
        new_start = end
        new_end = min(chr_len, end + upstream)
    else:
        print(f"Warning: unknown strand '{strand}', skipping.", file=sys.stderr)
        return None

    if new_start >= new_end:
        return None

    return chrom, new_start, new_end

def main():
    parser = argparse.ArgumentParser(description='Generate upstream 2kb BED from input BED and .fai.')
    parser.add_argument('bed', help='Input BED file (6 columns, with strand in 6th col)')
    parser.add_argument('fai', help='Genome index file (.fai)')
    parser.add_argument('-o', '--output', required=True, help='Output BED file')
    parser.add_argument('-u', '--upstream', type=int, default=2000,
                        help='Upstream distance in bp (default: 2000)')
    args = parser.parse_args()

    chrom_lengths = read_fai(args.fai)

    with open(args.bed, 'r') as fin, open(args.output, 'w') as fout:
        for line_num, line in enumerate(fin, 1):
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('track'):
                continue
            fields = line.split()
            if len(fields) < 6:
                print(f"Warning: line {line_num} has less than 6 fields, skipping.", file=sys.stderr)
                continue
            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                print(f"Warning: invalid start/end at line {line_num}, skipping.", file=sys.stderr)
                continue
            name = fields[3]
            score = fields[4]
            strand = fields[5]

            region = get_upstream_region(chrom, start, end, strand, chrom_lengths, args.upstream)
            if region is None:
                continue

            new_chrom, new_start, new_end = region
            fout.write(f"{new_chrom}\t{new_start}\t{new_end}\t{name}\t{score}\t{strand}\n")

    print(f"Done. Upstream BED written to {args.output}")

if __name__ == '__main__':
    main()
