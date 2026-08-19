#!/usr/bin/env python3
"""
For each query in the alignment score file, keep only the alignment with the
highest total score, and output query, target family, and score.
"""

import sys
import argparse

def extract_family(name):
    """Extract the family name from a miRNA isoform ID (e.g., MIR1030-isoform2-7 -> MIR1030)."""
    return name.split('-')[0]

def main():
    parser = argparse.ArgumentParser(
        description='Keep the best scoring alignment per query.')
    parser.add_argument('-i', '--input', required=True, help='Input score file')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    args = parser.parse_args()

    best = {}   # query -> (family, score)

    with open(args.input, 'r') as fin:
        first_line = fin.readline().strip()
        # If the first line looks like a header (starts with "query"), skip it
        if first_line.startswith("query") or first_line.startswith("#"):
            pass   # header skipped
        else:
            # Treat first line as data
            parts = first_line.split('\t')
            if len(parts) >= 6:
                query = parts[0]
                reference = parts[1]
                total_score = float(parts[5])
                family = extract_family(reference)
                if query not in best or total_score > best[query][1]:
                    best[query] = (family, total_score)

        # Process remaining lines
        for line in fin:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 6:
                continue
            query = parts[0]
            reference = parts[1]
            total_score = float(parts[5])
            family = extract_family(reference)
            if query not in best or total_score > best[query][1]:
                best[query] = (family, total_score)

    with open(args.output, 'w') as fout:
        for query in sorted(best.keys()):
            family, score = best[query]
            fout.write(f"{query}\t{family}\t{score:.1f}\n")

if __name__ == '__main__':
    main()
