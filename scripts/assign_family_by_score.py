#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
assign_family_by_score.py

Assign a single miRNA family to each query based on its best alignment score,
with special rules for certain family pairs (MIR156/MIR157, MIR165/MIR166,
MIR170/MIR171) that ignore scores and use sequence features instead.

Input format (tab-separated, 3 columns, no header):
    query    family    score

Output:
    -o : main output with chosen family (two columns: query, chosen_family)
    -s : secondary output containing all lines for queries that had multiple
         distinct families in the input.

Usage:
    python assign_family_by_score.py -i input.txt -o chosen.tsv -s multi.tsv [-f sequences.fasta]
"""

import sys
import argparse
import re
from collections import defaultdict

def extract_number(family):
    """Extract numeric part from family string, e.g., 'MIR1030' -> 1030."""
    match = re.search(r'\d+', family)
    if match:
        return int(match.group())
    return float('inf')

def read_fasta(fasta_file):
    """
    Read FASTA file and return a dictionary {header: sequence}.
    Header is taken as the entire line after '>' (without leading/trailing spaces).
    Sequence lines are concatenated.
    """
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
                    if current_header is not None:
                        seq_dict[current_header] = ''.join(current_seq)
                    current_header = line[1:].strip()
                    current_seq = []
                else:
                    current_seq.append(line.upper())
            if current_header is not None:
                seq_dict[current_header] = ''.join(current_seq)
    except IOError as e:
        sys.stderr.write(f"Error reading FASTA file: {e}\n")
        sys.exit(1)
    return seq_dict

def resolve_special_pair(f1, f2, query_seq, families_present, family_scores):
    """
    Attempt to use special rules for a conflicting pair (f1, f2).
    Returns the chosen family (f1 or f2) if the rule applies, otherwise None.
    family_scores is the {family: best_score} dict for the current query.
    """
    # MIR165 vs MIR166
    if {f1, f2} == {'MIR165', 'MIR166'} and query_seq is not None:
        if len(query_seq) >= 17:
            base = query_seq[16].upper()   # 17th base -> index 16
            if base == 'C':
                return 'MIR165'
            elif base == 'T':
                return 'MIR166'
        return None  # fallback if conditions not met
    # MIR170 vs MIR171
    if {f1, f2} == {'MIR170', 'MIR171'} and query_seq is not None:
        if len(query_seq) >= 12:
            base = query_seq[11].upper()   # 12th base -> index 11
            if base == 'C':
                return 'MIR170'
            elif base == 'T':
                return 'MIR171'
        return None
    # MIR156 vs MIR157
    if {f1, f2} == {'MIR156', 'MIR157'} and query_seq is not None:
        seq_len = len(query_seq)
        if seq_len == 20:
            if seq_len >= 11:
                base = query_seq[10].upper()   # 11th base -> index 10
                if base == 'G':
                    return 'MIR156'
                elif base == 'T':
                    return 'MIR157'
        elif seq_len == 21:
            if seq_len >= 12:
                base = query_seq[11].upper()   # 12th base -> index 11
                if base == 'G':
                    return 'MIR156'
                elif base == 'T':
                    return 'MIR157'
        return None
    return None

def decide_family(family_scores, query_seq=None):
    """
    Decide the best family for a query given a dict of family -> best_score.
    Special rules are applied to resolve ties/conflicts for known pairs.
    """
    if not family_scores:
        return None

    # 1. Apply special rules to merge predefined pairs if both present.
    #    We repeatedly apply rules as long as a pair exists.
    #    After merging, we keep the higher score as the score for the chosen family.
    applied = True
    while applied:
        applied = False
        families = list(family_scores.keys())
        for f1 in families:
            for f2 in families:
                if f1 >= f2:
                    continue
                # Check if this pair is a special pair and both are present
                winner = resolve_special_pair(f1, f2, query_seq, set(families), family_scores)
                if winner is not None:
                    # Merge: winner keeps the max of the two scores
                    merged_score = max(family_scores[f1], family_scores[f2])
                    # Remove both, add winner
                    del family_scores[f1]
                    del family_scores[f2]
                    family_scores[winner] = merged_score
                    applied = True
                    break
            if applied:
                break

    # 2. Find the maximum score among the remaining families
    max_score = max(family_scores.values())
    best_families = [fam for fam, scr in family_scores.items() if scr == max_score]

    # 3. If only one, use it. If multiple, choose the one with the smallest numeric part.
    best_families.sort(key=extract_number)
    return best_families[0]

def main():
    parser = argparse.ArgumentParser(
        description='Assign a single miRNA family to each query based on alignment scores, '
                    'with special sequence‑based rules for MIR156/157, MIR165/166 and MIR170/171.')
    parser.add_argument('-i', '--input', required=True,
                        help='Input file (query, family, score) without header')
    parser.add_argument('-o', '--output', required=True,
                        help='Output file with chosen family (two columns)')
    parser.add_argument('-s', '--secondary', required=True,
                        help='Output file for multi‑family queries (all original lines)')
    parser.add_argument('-f', '--fasta', default=None,
                        help='FASTA file with query sequences (required for special rules)')
    args = parser.parse_args()

    # Read all lines and group by query
    query_data = defaultdict(list)   # query -> list of (family, score, full_line)
    with open(args.input, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            query = parts[0].strip()
            family = parts[1].strip()
            try:
                score = float(parts[2])
            except ValueError:
                continue
            query_data[query].append((family, score, line))

    # Load sequences if FASTA provided
    seq_db = {}
    if args.fasta:
        seq_db = read_fasta(args.fasta)

    # Prepare output files
    fout_main = open(args.output, 'w')
    fout_sec = open(args.secondary, 'w')

    for query in sorted(query_data.keys()):
        entries = query_data[query]

        # Build family -> best score
        family_scores = {}
        for fam, scr, _ in entries:
            if fam not in family_scores or scr > family_scores[fam]:
                family_scores[fam] = scr

        # Decide best family
        chosen_family = decide_family(family_scores, seq_db.get(query))
        if chosen_family is None:
            chosen_family = "Unknown"

        fout_main.write(f"{query}\t{chosen_family}\n")

        # Write to secondary file if there are at least two distinct families
        distinct_fams = set(fam for fam, _, _ in entries)
        if len(distinct_fams) > 1:
            for _, _, line in entries:
                fout_sec.write(line + '\n')

    fout_main.close()
    fout_sec.close()

if __name__ == "__main__":
    main()
