#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
Split miRNA FASTA sequences into host and guest sets based on family size.

For each family (numeric ID >= -s/--start):
  - <=5 members : 1 host (rest guest)
  - 6-10 members: 2 hosts (rest guest)
  - >10 members : 3 hosts (rest guest)
Selection of hosts: if a -m mapping file is provided, prefer higher weight
(third column); ties are broken by larger count (suffix number), then longer
sequence, then random.  If a sequence has no weight, it is treated as 0.

Families with ID < start go entirely to guest.
"""

import argparse
import sys
import re
from collections import defaultdict
import random

def parse_fasta(filepath):
    """Read FASTA and return list of (header, sequence)."""
    seqs = []
    current_header = None
    current_seq = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header is not None:
                    seqs.append((current_header, ''.join(current_seq)))
                current_header = line[1:]   # remove '>'
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_header is not None:
            seqs.append((current_header, ''.join(current_seq)))
    return seqs

def parse_header(header):
    """
    Extract family and count from header like 'MIR12336-rename1-1'.
    Returns (family, count) or None if format doesn't match.
    """
    m = re.match(r'^(MIRN?\d+)-rename\d+-(\d+)$', header)
    if m:
        return m.group(1), int(m.group(2))
    return None

def family_numeric_id(family_str):
    """Extract integer from family string like 'MIR12336'."""
    m = re.search(r'\d+', family_str)
    return int(m.group()) if m else 0

def load_weights(weight_file):
    """
    Read the -m file (two or three columns, tab-separated).
    Returns a dict: header -> weight (float). Missing third column -> 0.0.
    """
    weights = {}
    with open(weight_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            hdr = parts[0].strip()
            if len(parts) >= 3 and parts[2].strip():
                try:
                    w = float(parts[2])
                except ValueError:
                    sys.stderr.write(f"Warning: invalid weight '{parts[2]}' for {hdr}, set to 0.\n")
                    w = 0.0
            else:
                w = 0.0
            weights[hdr] = w
    return weights

def select_hosts(members, num_hosts):
    """
    members: list of (header, seq, count, weight)
    Select exactly num_hosts hosts (or all if less).
    Sorting: weight desc, count desc, length desc, then random.
    """
    n = len(members)
    take = min(num_hosts, n)
    if take == 0:
        return []
    # Sort by weight desc, count desc, length desc
    members_sorted = sorted(members, key=lambda x: (-x[3], -x[2], -len(x[1])))
    # Break ties randomly within groups that have identical weight, count, and length
    selected = []
    i = 0
    while i < len(members_sorted) and len(selected) < take:
        j = i
        while j < len(members_sorted) and \
              members_sorted[j][3] == members_sorted[i][3] and \
              members_sorted[j][2] == members_sorted[i][2] and \
              len(members_sorted[j][1]) == len(members_sorted[i][1]):
            j += 1
        group = members_sorted[i:j]
        random.shuffle(group)
        for item in group:
            if len(selected) < take:
                selected.append(item)
            else:
                break
        i = j
    return selected

def main():
    parser = argparse.ArgumentParser(
        description='Select host sequences per family based on count, length, and optional weight.')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-s', '--start', type=int, required=True,
                        help='Numeric threshold: families with ID >= start are processed.')
    parser.add_argument('-m', '--weight', default=None,
                        help='Optional weight file (two or three columns: name, family, weight)')
    parser.add_argument('--host', required=True, help='Output host FASTA')
    parser.add_argument('--guest', required=True, help='Output guest FASTA')
    args = parser.parse_args()

    random.seed(42)  # reproducibility

    sequences = parse_fasta(args.input)

    # Load weights if provided
    weights = {}
    if args.weight:
        weights = load_weights(args.weight)

    # Group by family, storing also the weight for each member
    family_map = defaultdict(list)  # family -> [(header, seq, count, weight), ...]
    unparsed = []  # headers that didn't match pattern -> guest

    for hdr, seq in sequences:
        info = parse_header(hdr)
        if info:
            family, count = info
            w = weights.get(hdr, 0.0)   # default 0 if not in weight file
            family_map[family].append((hdr, seq, count, w))
        else:
            unparsed.append((hdr, seq))

    host_entries = []
    guest_entries = []

    for family, members in family_map.items():
        fid = family_numeric_id(family)
        # If below threshold, entire family goes to guest
        if fid < args.start:
            for hdr, seq, _, _ in members:
                guest_entries.append((hdr, seq))
            continue

        n = len(members)
        if n <= 5:
            num_hosts = 1
        elif n <= 10:
            num_hosts = 2
        else:
            num_hosts = 3

        selected = select_hosts(members, num_hosts)
        selected_headers = {hdr for hdr, _, _, _ in selected}

        for hdr, seq, _, _ in members:
            if hdr in selected_headers:
                host_entries.append((hdr, seq))
            else:
                guest_entries.append((hdr, seq))

    # Unparsed headers go to guest
    for hdr, seq in unparsed:
        guest_entries.append((hdr, seq))

    # Write outputs
    with open(args.host, 'w') as fh:
        for hdr, seq in host_entries:
            fh.write(f">{hdr}\n{seq}\n")

    with open(args.guest, 'w') as fg:
        for hdr, seq in guest_entries:
            fg.write(f">{hdr}\n{seq}\n")

    print(f"Host sequences: {len(host_entries)}", file=sys.stderr)
    print(f"Guest sequences: {len(guest_entries)}", file=sys.stderr)

if __name__ == '__main__':
    main()
