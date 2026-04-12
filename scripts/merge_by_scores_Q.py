#!/usr/bin/env python3
"""
merge_by_query.py

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

Usage:
    python merge_by_query.py -i isoform.best.score -o merged.txt
    cat isoform.best.score | python merge_by_query.py > merged.txt
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

def select_best_family(family_counts):
    """
    family_counts: dict {family: count}
    Select best family based on:
        1. Highest count (descending)
        2. Tie: family_priority_key (ascending)
        3. If still tie: alphabetical order
    """
    best_family = max(family_counts.items(), key=lambda x: (x[1], -family_priority_key(x[0])[0], -family_priority_key(x[0])[1] if family_priority_key(x[0])[1] != float('inf') else 0))
    # 上面 max 的 key 需要仔细处理：先按 count 降序，再按 priority 升序，再按 number 升序。
    # 由于 max 默认是升序，我们使用负数来达到降序效果。更清晰的方式是用 sorted：
    # 这里改用 sorted 以便理解：
    # 但为了简洁，重新实现：
    # 实际上直接写清楚：
    items = list(family_counts.items())
    # 排序：count 降序，priority 升序，number 升序，family 字母升序
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

def main():
    parser = argparse.ArgumentParser(description='Merge best score file by query, summarize target families per query.')
    parser.add_argument('-i', '--input', help='Input best score file (default: stdin)')
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
        for family in sorted(family_counts.keys()):  # alphabetical order
            count = family_counts[family]
            if count == 1:
                summary_parts.append(family)
            else:
                summary_parts.append(f"{family}-{count}")
        summary_str = ','.join(summary_parts)

        # Select best family for column 3
        selected_family = select_best_family(family_counts)

        outfile.write(f"{query}\t{summary_str}\t{selected_family}\n")

    if outfile is not sys.stdout:
        outfile.close()

if __name__ == "__main__":
    main()
