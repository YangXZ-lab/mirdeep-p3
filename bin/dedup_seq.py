#!/usr/bin/env python3
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
Global sequence-level deduplication of a FASTA file.

The first occurrence of each unique sequence is written to the unique file;
all subsequent duplicates are written to the duplicates file.
Matching is based on sequence identity only (case-insensitive), ignoring headers.

Usage:
    python3 dedup_seq.py -i input.fa --uni unique.fa --dup duplicates.fa
"""

import argparse
import sys


def read_fasta(path):
    """Read FASTA file, return list of (header, seq)."""
    records = []
    with open(path) as f:
        header = None
        seq_parts = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
        if header is not None:
            records.append((header, "".join(seq_parts)))
    return records


def main():
    parser = argparse.ArgumentParser(
        description="Global FASTA dedup: first occurrence -> unique, rest -> duplicates")
    parser.add_argument("-i", "--input", required=True,
                        help="Input FASTA file")
    parser.add_argument("--uni", required=True,
                        help="Output file for unique sequences")
    parser.add_argument("--dup", required=True,
                        help="Output file for duplicate sequences")
    args = parser.parse_args()

    sys.stderr.write(f"Reading: {args.input}\n")
    records = read_fasta(args.input)
    sys.stderr.write(f"  total records: {len(records)}\n")

    seen = set()
    uniq_lines = []
    dup_lines = []

    for hdr, seq in records:
        if seq in seen:
            dup_lines.append(f">{hdr}\n{seq}\n")
        else:
            seen.add(seq)
            uniq_lines.append(f">{hdr}\n{seq}\n")

    with open(args.uni, "w") as f:
        f.writelines(uniq_lines)
    with open(args.dup, "w") as f:
        f.writelines(dup_lines)

    n_total = len(records)
    n_uniq = len(uniq_lines)
    n_dup = len(dup_lines)
    pct = n_uniq / n_total * 100 if n_total else 0

    sys.stderr.write(f"\n===== Summary =====\n")
    sys.stderr.write(f"Total records:  {n_total}\n")
    sys.stderr.write(f"Unique:         {n_uniq} ({pct:.1f}%)\n")
    sys.stderr.write(f"Duplicates:     {n_dup} ({100-pct:.1f}%)\n")
    sys.stderr.write(f"Unique file:    {args.uni}\n")
    sys.stderr.write(f"Duplicates file: {args.dup}\n")


if __name__ == "__main__":
    main()
