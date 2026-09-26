#!/usr/bin/env python3
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
Deduplicate alignment score rows by unordered (query, reference) pair.

Rows sharing the same unordered pair of column 1 and column 2 are grouped;
only rows with the maximum total_score (column 6) are kept. Ties are all
retained (e.g. both forward and reverse rows with the same top score).

Usage:
    python3 dedup_score_pair.py -i input.tsv -o output.tsv
"""

import argparse
import sys
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(
        description="Keep max-score row per unordered pair (col1, col2)")
    parser.add_argument("-i", "--input", required=True,
                        help="Input TSV file")
    parser.add_argument("-o", "--output", required=True,
                        help="Output TSV file")
    args = parser.parse_args()

    groups = defaultdict(list)
    n_total = 0
    header = None

    with open(args.input) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if header is None and not line[0].isdigit():
                header = line
                continue
            cols = line.split("\t")
            if len(cols) < 6:
                continue
            q, r = cols[0], cols[1]
            key = (q, r) if q <= r else (r, q)
            groups[key].append((float(cols[5]), line))
            n_total += 1

    n_out = 0
    with open(args.output, "w") as f:
        if header is not None:
            f.write(header + "\n")
        for key in sorted(groups.keys()):
            rows = groups[key]
            best = max(score for score, _ in rows)
            for score, line in rows:
                if score == best:
                    f.write(line + "\n")
                    n_out += 1

    sys.stderr.write(f"Total rows:      {n_total}\n")
    sys.stderr.write(f"Unique pairs:    {len(groups)}\n")
    sys.stderr.write(f"Output rows:     {n_out} (best-score rows)\n")
    sys.stderr.write(f"Removed:         {n_total - n_out}\n")
    sys.stderr.write(f"Output: {args.output}\n")


if __name__ == "__main__":
    main()
