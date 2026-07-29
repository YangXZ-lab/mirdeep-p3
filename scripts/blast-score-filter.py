#!/usr/bin/env python3
"""
Filter alignment scores based on per‑query maximum and user thresholds.

Usage:
    python filter_scores.py -s scores.tsv -m max.tsv -t 55 -d 20 -o out.tsv
"""

import argparse
import sys
from collections import defaultdict

def extract_family(name):
    """Extract the family part from a miRNA isoform name (e.g., 'MIR472' from 'MIR472-isoform16-5')."""
    return name.split('-')[0]

def main():
    parser = argparse.ArgumentParser(description='Filter score file by query maximum and thresholds.')
    parser.add_argument('-s', '--scores', required=True,
                        help='Input score file (6 columns, with header)')
    parser.add_argument('-m', '--max', required=True,
                        help='Per‑query max score file (3 columns, with or without header)')
    parser.add_argument('-t', '--threshold', type=float, default=55.0,
                        help='Absolute score threshold (default: 55)')
    parser.add_argument('-d', '--delta', type=float, default=20.0,
                        help='Allowed drop from maximum (default: 20)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output filtered file')
    args = parser.parse_args()

    # 1. Read per-query maximum scores (skip header if present)
    qmax = {}
    with open(args.max, 'r') as f:
        first = True
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Guess if header by looking at first line
            if first:
                first = False
                # If the first line starts with typical header words, skip
                if line.startswith('query') or 'max' in line.lower():
                    continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            query = parts[0].strip()
            try:
                score = float(parts[2])
            except ValueError:
                continue
            qmax[query] = score

    # Keep only queries whose max score >= threshold
    valid_queries = {q for q, s in qmax.items() if s >= args.threshold}

    # 2. Aggregate best score per (query, family)
    best = defaultdict(float)   # (query, family) -> max_score
    with open(args.scores, 'r') as f:
        first = True
        for line in f:
            line = line.strip()
            if not line:
                continue
            if first:
                first = False
                if line.startswith('query') or 'total_score' in line:
                    continue
            parts = line.split('\t')
            if len(parts) < 6:
                continue
            query = parts[0].strip()
            if query not in valid_queries:
                continue
            reference = parts[1].strip()
            try:
                total = float(parts[5])
            except ValueError:
                continue
            family = extract_family(reference)
            key = (query, family)
            if total > best[key]:
                best[key] = total

    # 3. Apply threshold and dynamic drop
    with open(args.output, 'w') as out:
        # Sort output by query then family for reproducibility
        for (query, family), score in sorted(best.items(), key=lambda x: (x[0][0], x[0][1])):
            max_q = qmax[query]
            cutoff = max(args.threshold, max_q - args.delta)
            if score >= cutoff:
                out.write(f"{query}\t{family}\t{score:.1f}\n")

if __name__ == '__main__':
    main()
