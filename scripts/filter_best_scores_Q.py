#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
filter_best_scores_Q.py

Filter bowtie scoring results to keep only the best scoring alignment per query.
Best is defined as: minimum total score (S). If total scores tie, then compare
subscores in order: B, A, C, D (lower is better). If still tied, compare
family priority (MIR < MIRN < other) and family number (smaller better).
If all are identical, keep all such lines.

Family extraction: from query name before '-isoform' (if present), else whole name.
MIR family: starts with 'MIR' but not 'MIRN'.
MIRN family: starts with 'MIRN'.
Others: any other pattern.

Self-alignments are not excluded a priori; they compete under the same rules.

Input: tab-separated with columns:
    query target strand A B C D total
Output: same format, filtered.

Usage:
    python filter_best_scores_Q.py -i scores.txt -o best_scores.txt
    cat scores.txt | python filter_best_scores_Q.py > best_scores.txt
"""

import sys
import argparse
import re
from collections import defaultdict

def extract_family(query):
    """Extract miRNA family from query name (before '-isoform')."""
    if '-isoform' in query:
        return query.split('-isoform')[0]
    return query

def get_family_priority_and_number(query):
    """
    Return (priority, number) for sorting.
    priority: 0 for MIR, 1 for MIRN, 2 for others.
    number: extracted integer from family, or inf if none.
    """
    family = extract_family(query)
    # Determine type
    if family.startswith('MIR') and not family.startswith('MIRN'):
        priority = 0
    elif family.startswith('MIRN'):
        priority = 1
    else:
        priority = 2

    # Extract number
    match = re.search(r'\d+', family)
    if match:
        number = int(match.group())
    else:
        number = float('inf')
    return (priority, number)

def read_scores(file_handle):
    """Generator yielding tuples (query, target, strand, A, B, C, D, total) from file."""
    for line in file_handle:
        line = line.strip()
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) != 8:
            sys.stderr.write(f"Warning: skipping invalid line (not 8 columns): {line}\n")
            continue
        query, target, strand, A, B, C, D, total = parts
        try:
            A = float(A)
            B = float(B)
            C = float(C)
            D = float(D)
            total = float(total)
        except ValueError:
            sys.stderr.write(f"Warning: skipping line with non-numeric scores: {line}\n")
            continue
        yield (query, target, strand, A, B, C, D, total)

def score_key(entry):
    """
    Return a tuple for sorting:
    (total, B, A, C, D, family_priority, family_number) all ascending.
    """
    query, _, _, A, B, C, D, total = entry
    priority, number = get_family_priority_and_number(query)
    return (total, B, A, C, D, priority, number)

def filter_best_per_query(entries):
    """
    Given an iterable of entries, group by query and yield the best ones
    according to the scoring rules.
    """
    groups = defaultdict(list)
    for e in entries:
        groups[e[0]].append(e)   # group by query

    for query, group in groups.items():
        sorted_group = sorted(group, key=score_key)
        best_key = score_key(sorted_group[0])
        for entry in sorted_group:
            if score_key(entry) == best_key:
                yield entry
            else:
                break

def main():
    parser = argparse.ArgumentParser(description='Filter bowtie scoring results to keep best per query with family priority.')
    parser.add_argument('-i', '--input', help='Input scores file (default: stdin)')
    parser.add_argument('-o', '--output', help='Output file (default: stdout)')
    args = parser.parse_args()

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

    # Read entries
    entries = list(read_scores(infile))
    if infile is not sys.stdin:
        infile.close()

    # Filter
    best_entries = list(filter_best_per_query(entries))

    # Write output
    for e in best_entries:
        query, target, strand, A, B, C, D, total = e
        outfile.write(f"{query}\t{target}\t{strand}\t{A:.1f}\t{B:.1f}\t{C:.1f}\t{D:.1f}\t{total:.1f}\n")

    if outfile is not sys.stdout:
        outfile.close()

if __name__ == "__main__":
    main()
