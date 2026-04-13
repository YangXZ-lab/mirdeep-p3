#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.


"""
merge_by_scores_Q.py

Merge the best score file by query, collecting unique target families and their counts.
Input format (tab-separated, 8 columns):
    query target strand A B C D total
Output format (tab-separated, 3 columns):
    query    family_summary    selected_family

Family summary format: for each family, if count=1 -> "family", else "family-count".
Families are sorted alphabetically. Example: "MIR1030-2,MIR1130"

Selection rule for third column (based on families):
    1. Choose family with highest count.
    2. If tie, prefer MIR > MIRN > other, and smaller numeric part.
    3. If still tie, alphabetically first.

SPECIAL RULE (added):
    If the candidate families contain both 'MIR165' and 'MIR166',
    then the decision is made based on the 17th base of the query sequence:
        - 17th base == 'C' -> choose MIR165
        - 17th base == 'T' -> choose MIR166
    This requires providing a FASTA file via -f/--fasta.

Usage:
    python merge_by_query.py -i isoform.best.score -o merged.txt -f sequences.fasta
    cat isoform.best.score | python merge_by_query.py -f sequences.fasta > merged.txt
"""

import sys
import argparse
import re
from collections import defaultdict

def extract_family(name):
    """
    Extract family from a name (query or target) by taking part before '-isoform'.
    Example: 'MIR1030-isoform4-1' -> 'MIR1030'
             'MIR1130-isoform2-1' -> 'MIR1130'
    """
    if '-isoform' in name:
        return name.split('-isoform')[0]
    return name

def extract_number(family):
    """Extract numeric part from family string, e.g., 'MIR1030' -> 1030."""
    match = re.search(r'\d+', family)
    if match:
        return int(match.group())
    return float('inf')

def family_priority_key(family):
    """
    Return (priority, number) for sorting families.
    priority: 0 for MIR (starts with 'MIR' but not 'MIRN'),
              1 for MIRN,
              2 for others.
    number: extracted integer (lower is better).
    """
    if family.startswith('MIR') and not family.startswith('MIRN'):
        priority = 0
    elif family.startswith('MIRN'):
        priority = 1
    else:
        priority = 2
    return (priority, extract_number(family))

def select_best_family(family_counts, query_seq=None):
    """
    family_counts: dict {family: count}
    query_seq: optional string, the sequence of the query (for special 165/166 rule)
    
    Select best family based on:
        1. SPECIAL: if both 'MIR165' and 'MIR166' in family_counts and query_seq is provided,
           use 17th base: C -> MIR165, T -> MIR166.
        2. Otherwise: highest count -> priority -> number -> alphabetical.
    """
    families = set(family_counts.keys())
    
    # Special handling for MIR165 vs MIR166 conflict
    if 'MIR165' in families and 'MIR166' in families and query_seq is not None:
        if len(query_seq) >= 17:
            base = query_seq[16].upper()  # 17th base (1-indexed) -> index 16
            if base == 'C':
                return 'MIR165'
            elif base == 'T':
                return 'MIR166'
            else:
                sys.stderr.write(f"Warning: 17th base is '{base}', not C or T. Falling back to default rule.\n")
        else:
            sys.stderr.write(f"Warning: Query sequence length {len(query_seq)} < 17. Falling back to default rule.\n")
    
    # Default rule
    items = list(family_counts.items())
    # Sort: count descending, priority ascending, number ascending, family alphabetical
    items.sort(key=lambda x: (-x[1], family_priority_key(x[0])[0], family_priority_key(x[0])[1], x[0]))
    return items[0][0]

def read_best_file(file_handle):
    """Yield (query, target) for each line in file."""
    for line in file_handle:
        line = line.strip()
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) != 8:
            sys.stderr.write(f"Warning: skipping invalid line (not 8 columns): {line}\n")
            continue
        query, target, _, _, _, _, _, _ = parts
        yield query, target

def read_fasta(fasta_file):
    seq_dict = {}
    current_header = None
    current_seq = []
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    # Save previous sequence
                    if current_header is not None:
                        seq_dict[current_header] = ''.join(current_seq)
                    current_header = line[1:].strip()  # remove '>' and trim spaces
                    current_seq = []
                else:
                    current_seq.append(line)
            # Don't forget the last sequence
            if current_header is not None:
                seq_dict[current_header] = ''.join(current_seq)
    except IOError as e:
        sys.stderr.write(f"Error reading FASTA file: {e}\n")
        sys.exit(1)
    return seq_dict

def main():
    parser = argparse.ArgumentParser(description='Merge best score file by query, summarize target families per query.')
    parser.add_argument('-i', '--input', help='Input best score file (default: stdin)')
    parser.add_argument('-o', '--output', help='Output file (default: stdout)')
    parser.add_argument('-f', '--fasta', help='FASTA file containing query sequences (required for special 165/166 rule)')
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

    # Load FASTA if provided
    seq_db = {}
    if args.fasta:
        seq_db = read_fasta(args.fasta)
        sys.stderr.write(f"Loaded {len(seq_db)} sequences from {args.fasta}\n")

    # Group targets by query
    query_to_targets = defaultdict(list)
    for query, target in read_best_file(infile):
        query_to_targets[query].append(target)

    if infile is not sys.stdin:
        infile.close()

    # Process each query
    for query in sorted(query_to_targets.keys()):
        targets = query_to_targets[query]
        # Extract families and count
        family_counts = defaultdict(int)
        for t in targets:
            family = extract_family(t)
            family_counts[family] += 1

        # Build summary string for column 2
        summary_parts = []
        for family in sorted(family_counts.keys()):
            count = family_counts[family]
            if count == 1:
                summary_parts.append(family)
            else:
                summary_parts.append(f"{family}-{count}")
        summary_str = ','.join(summary_parts)

        # Retrieve query sequence for special rule (if available)
        query_seq = seq_db.get(query)  # None if not found or no FASTA provided
        
        # Select best family (pass sequence if we have it)
        selected_family = select_best_family(family_counts, query_seq)

        outfile.write(f"{query}\t{summary_str}\t{selected_family}\n")

    if outfile is not sys.stdout:
        outfile.close()

if __name__ == "__main__":
    main()
