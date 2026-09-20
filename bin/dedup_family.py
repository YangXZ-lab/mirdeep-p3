#!/usr/bin/env python3
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
Deduplicate FASTA sequences by miRNA family.

Families are identified by their numeric portion only. Both 'miR156a' and
'miR156d*' style headers (trailing letters, '*', or backtick) are merged into
'miR156'. Within each family, identical sequences are merged and labelled as
isoform1, isoform2, ... with a count suffix indicating how many copies existed
before dedup.

Usage:
    python3 dedup_family.py -i input.fa -o output.fa
"""

import argparse
import re
import sys
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(
        description="Deduplicate FASTA by miRNA family with isoform numbering")
    parser.add_argument("-i", "--input", required=True,
                        help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True,
                        help="Output FASTA file")
    args = parser.parse_args()

    # Group sequences by family
    family_seqs = defaultdict(list)
    family = None
    seq_parts = []

    with open(args.input) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if family is not None and seq_parts:
                    family_seqs[family].append("".join(seq_parts))
                # Extract family number: "Ach-miR156a" -> "miR156",
                # "Cre-MIRN10192a" -> "MIRN10192", "Smo-miR156d*" -> "miR156"
                m = re.search(r"(mi[Rr]N?\d+)", line)
                family = m.group(1) if m else re.sub(r"[a-z]+$", "", line[1:])
                seq_parts = []
            else:
                seq_parts.append(line)
        if family is not None and seq_parts:
            family_seqs[family].append("".join(seq_parts))

    total_in = sum(len(v) for v in family_seqs.values())
    total_out = 0

    with open(args.output, "w") as f:
        for fam in sorted(family_seqs.keys()):
            cnt = defaultdict(int)
            for s in family_seqs[fam]:
                cnt[s] += 1
            sorted_items = sorted(cnt.items(), key=lambda x: (-x[1], x[0]))
            for i, (seq, c) in enumerate(sorted_items, 1):
                f.write(f">{fam}-isoform{i}-{c}\n{seq}\n")
                total_out += 1

    sys.stderr.write(f"Input records:  {total_in}\n")
    sys.stderr.write(f"Output isoforms: {total_out} ({total_in - total_out} removed)\n")
    sys.stderr.write(f"Families:       {len(family_seqs)}\n")
    sys.stderr.write(f"Output file:    {args.output}\n")


if __name__ == "__main__":
    main()
