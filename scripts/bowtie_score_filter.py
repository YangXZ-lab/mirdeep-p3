#!/usr/bin/env python3
"""
Filter pairwise alignment scores: for each query, keep rows with the highest
total score. In case of ties, prefer lower values of score1, then score2,
then score3. If all scores are identical, all tied rows are kept.
"""

import sys
import argparse
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(
        description='Keep best scoring reference(s) per query based on total and tie-breaking rules.')
    parser.add_argument('-i', '--input', required=True, help='Input file (tab-separated).')
    parser.add_argument('-o', '--output', required=True, help='Output filtered file.')
    args = parser.parse_args()

    # Read all lines and group by query
    groups = defaultdict(list)
    with open(args.input, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 6:
                sys.stderr.write(f"Warning: skipping line with <6 columns: {line}\n")
                continue
            query = parts[0]
            try:
                s1 = float(parts[2])
                s2 = float(parts[3])
                s3 = float(parts[4])
                total = float(parts[5])
            except ValueError:
                sys.stderr.write(f"Warning: non-numeric score in line: {line}\n")
                continue
            groups[query].append((line, s1, s2, s3, total))

    # Process each group
    with open(args.output, 'w') as out:
        for query, rows in groups.items():
            if not rows:
                continue
            # Find the maximum total
            max_total = max(row[4] for row in rows)
            # Keep only those with the max total
            candidates = [r for r in rows if r[4] == max_total]
            # Tie-breaking: keep rows with minimum score1
            min_s1 = min(r[1] for r in candidates)
            candidates = [r for r in candidates if r[1] == min_s1]
            # Minimum score2
            min_s2 = min(r[2] for r in candidates)
            candidates = [r for r in candidates if r[2] == min_s2]
            # Minimum score3
            min_s3 = min(r[3] for r in candidates)
            candidates = [r for r in candidates if r[3] == min_s3]
            # Write all remaining (ties in all scores are kept)
            for r in candidates:
                out.write(r[0] + '\n')

if __name__ == '__main__':
    main()
