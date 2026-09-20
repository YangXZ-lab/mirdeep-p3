#!/usr/bin/env python3
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
Remove sequences from -i that have exact matches in -r (same length & sequence).

Builds a reference set from -r, then writes -i records whose sequence is NOT
present in the reference. Matching is case-insensitive and compares full-length
sequence identity.

Usage:
    python3 remove_matched_seq.py -i input.fa -r reference.fa -o output.fa
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
        description="Remove exact-match sequences from -i that appear in -r")
    parser.add_argument("-i", "--input", required=True,
                        help="Input FASTA file")
    parser.add_argument("-r", "--reference", required=True,
                        help="Reference FASTA file (sequences to remove)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output FASTA file")
    args = parser.parse_args()

    # Build reference set
    sys.stderr.write(f"Reading reference: {args.reference}\n")
    ref_records = read_fasta(args.reference)
    ref_seqs = {seq for _, seq in ref_records}
    sys.stderr.write(f"  unique ref sequences: {len(ref_seqs)}\n")

    # Filter input
    sys.stderr.write(f"Reading input: {args.input}\n")
    input_records = read_fasta(args.input)

    kept = 0
    removed = 0
    with open(args.output, "w") as out:
        for hdr, seq in input_records:
            if seq in ref_seqs:
                removed += 1
            else:
                out.write(f">{hdr}\n{seq}\n")
                kept += 1

    n_total = len(input_records)
    sys.stderr.write(f"Total input:  {n_total}\n")
    sys.stderr.write(f"Kept:         {kept}\n")
    sys.stderr.write(f"Removed:      {removed}\n")
    sys.stderr.write(f"Output: {args.output}\n")


if __name__ == "__main__":
    main()
