#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.


"""
bowtie_mismatch_scoring.py

Parse bowtie alignment output (custom format) and compute penalty scores:
A = sum of p for mismatches outside positions 10/11 (1-based)
B = sum of p * 2 for mismatches at positions 10 or 11
C = L - distance between two mismatches/indels (if exactly two events)
D = sum of (L - |pos - 9.5|) for each event
Total = A + B + C + D

If query length differs from target length (provided by a length file):
    - If quality string (col 6) consists of identical characters (e.g., all 'I'),
      assume indels occur at edges (ignored), but still process any mismatches from col 8.
    - Otherwise, treat each non-identical position in quality string as an indel event
      (p=1) and merge with mismatches from col 8, then compute scores normally.

Output: query_name target_name strand A B C D total (tab-separated)

Usage:
    python bowtie_mismatch_scoring.py -i input.aln -l isoform.length -o output.txt
"""

import sys
import argparse

def read_length_file(length_file):
    """Read tab-separated file with isoform name and length, return dict {name: length}."""
    length_dict = {}
    try:
        with open(length_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    sys.stderr.write(f"Warning: skipping invalid length line: {line}\n")
                    continue
                name = parts[0]
                try:
                    length = int(parts[1])
                    length_dict[name] = length
                except ValueError:
                    sys.stderr.write(f"Warning: invalid length value in line: {line}\n")
    except IOError as e:
        sys.stderr.write(f"Error reading length file: {e}\n")
        sys.exit(1)
    return length_dict

def parse_mismatch_string(m_str):
    """Parse mismatch column like '8:A>G,15:G>A' into list of (pos, q_base, r_base)."""
    mismatches = []
    if not m_str:
        return mismatches
    parts = m_str.split(',')
    for part in parts:
        if not part:
            continue
        try:
            pos_str, change = part.split(':')
            q_base, r_base = change.split('>')
            pos = int(pos_str)
            mismatches.append((pos, q_base, r_base))
        except ValueError:
            sys.stderr.write(f"Warning: cannot parse mismatch '{part}'\n")
    return mismatches

def get_p(q_base, r_base):
    """Return p based on mismatch type: 0.5 for A<->G or C<->T, else 1.0."""
    if (q_base, r_base) in [('A','G'), ('G','A'), ('C','T'), ('T','C')]:
        return 0.5
    else:
        return 1.0

def all_same_char(s):
    """Return True if string s is empty or all characters are identical."""
    if not s:
        return True
    first = s[0]
    return all(ch == first for ch in s)

def process_line(line, length_dict, normal_qual_char='I'):
    """
    Process one alignment line.
    Returns tuple (query, target, strand, A, B, C, D) or None if line is invalid.
    """
    fields = line.strip().split('\t')
    # Need at least 6 columns: query, strand, target, _, seq, qual
    if len(fields) < 6:
        sys.stderr.write(f"Warning: line has too few columns, skipped: {line.strip()}\n")
        return None

    query = fields[0]
    strand = fields[1]
    target = fields[2]
    seq = fields[4]                # sequence
    qual = fields[5]                # quality string (used for indel detection)
    mism_str = fields[7] if len(fields) >= 8 else None

    # Get target length from dictionary
    if target not in length_dict:
        sys.stderr.write(f"Warning: target '{target}' not found in length file; assuming no indel.\n")
        target_len = len(seq)
    else:
        target_len = length_dict[target]

    query_len = len(seq)

    # Parse mismatches (if any)
    mismatches = parse_mismatch_string(mism_str) if mism_str else []
    mismatch_events = [(pos, get_p(qb, rb)) for pos, qb, rb in mismatches]

    # Determine events list
    if query_len == target_len:
        # No indel, only mismatches
        events = mismatch_events
    else:
        # Length differs, potential indels
        if all_same_char(qual):
            # Edge indels (ignored), only mismatches count
            events = mismatch_events
        else:
            # Middle indels: treat each non-normal position as indel event (p=1)
            indel_events = []
            for i, ch in enumerate(qual):
                if ch != normal_qual_char:
                    indel_events.append((i, 1.0))
            events = indel_events + mismatch_events
            # Optional warning
            if abs(query_len - target_len) != len(indel_events):
                sys.stderr.write(f"Warning: length difference ({abs(query_len - target_len)}) "
                                 f"does not match number of putative indels ({len(indel_events)}) "
                                 f"in line: {line.strip()}\n")

    # If no events, all scores zero
    if not events:
        return (query, target, strand, 0.0, 0.0, 0.0, 0.0)

    # Separate positions and p values
    positions = [pos for pos, _ in events]
    p_vals = [p for _, p in events]

    # Calculate A and B
    A = 0.0
    B = 0.0
    for pos, p in zip(positions, p_vals):
        if pos == 9 or pos == 10:   # 0-based positions for 10th and 11th bases
            B += p * 2
        else:
            A += p * 1

    # Calculate D
    center = 9.5
    D = 0.0
    for pos in positions:
        dist_from_center = abs(pos - center)
        D += (query_len - dist_from_center)

    # Calculate C if exactly two events
    C = 0.0
    if len(events) == 2:
        dist_between = abs(positions[0] - positions[1])
        C = query_len - dist_between
    elif len(events) > 2:
        sys.stderr.write(f"Warning: more than two events ({len(events)}) found; C set to 0. Line: {line.strip()}\n")

    return (query, target, strand, A, B, C, D)

def main():
    parser = argparse.ArgumentParser(description='Score bowtie alignments based on mismatch and indel penalties.')
    parser.add_argument('-i', '--input', help='Input alignment file (default: stdin)')
    parser.add_argument('-o', '--output', help='Output file (default: stdout)')
    parser.add_argument('-l', '--length', required=True, help='Length file (tab-separated: isoform_name length)')
    parser.add_argument('-n', '--normal-char', default='I', help='Character in quality string representing normal match (default: I)')
    args = parser.parse_args()

    # Read length file
    length_dict = read_length_file(args.length)

    # Open input
    if args.input:
        try:
            infile = open(args.input, 'r')
        except IOError as e:
            sys.stderr.write(f"Error opening input file: {e}\n")
            sys.exit(1)
    else:
        infile = sys.stdin

    # Open output
    if args.output:
        try:
            outfile = open(args.output, 'w')
        except IOError as e:
            sys.stderr.write(f"Error opening output file: {e}\n")
            sys.exit(1)
    else:
        outfile = sys.stdout

    for line in infile:
        line = line.rstrip('\n')
        if not line:
            continue
        result = process_line(line, length_dict, args.normal_char)
        if result:
            query, target, strand, A, B, C, D = result
            total = A + B + C + D
            outfile.write(f"{query}\t{target}\t{strand}\t{A:.1f}\t{B:.1f}\t{C:.1f}\t{D:.1f}\t{total:.1f}\n")

    if infile is not sys.stdin:
        infile.close()
    if outfile is not sys.stdout:
        outfile.close()

if __name__ == "__main__":
    main()
