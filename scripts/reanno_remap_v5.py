#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
Combine cluster, singleton and alignment results into a three-column
miRNA-to-family mapping file.

Third column (ratio) is preserved for entries from the --multi file;
for other sources it is left empty.

Naming behaviour:
  - If a cluster/singleton can reuse its intrinsic family number (derived
    from the miRNA name) and that family is not in the known FASTA set and
    not yet used, it is named with the fixed prefix "MIR".
  - If it cannot reuse (conflict or unknown), it is assigned a new sequential
    number using the prefix given by --type (MIR or MIRN).
  - --type MIR   : new families use "MIR" prefix; -s/--start is required.
  - --type MIRN  : new families use "MIRN" prefix; -s/--start must be omitted,
                   numbering starts at 1.

The input files --alncl, --single, --multi can be used in any combination:
  - all three provided,
  - only one or two of them provided.
At least one of them should be provided to generate a meaningful output.

Usage:
    python merge_family.py --alncl alncl.txt --single single.txt --multi multi.txt \
                           --type MIRN -f known.fasta -o output.txt
    python merge_family.py --single single.txt --type MIR -s 12100 -f known.fasta -o output.txt
    python merge_family.py --multi multi.txt --type MIRN -f known.fasta -o output.txt
"""

import argparse
import sys
import re
from collections import defaultdict, Counter

def parse_fasta_families(fasta_file):
    """
    Extract family numeric identifiers from FASTA headers.
    Supports both 'MIR' and 'MIRN' prefixes (e.g. >MIR2111... or >MIRN1...).
    Returns a set of numeric strings, e.g. {'2111', '1'}.
    """
    families = set()
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                m = re.match(r'>MIRN?(\d+)', line)
                if m:
                    families.add(m.group(1))
    return families

def parse_alncl(filepath):
    """Read alncl file (two columns: mirna, family) -> list of (mirna, family, '')"""
    entries = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                sys.stderr.write(f"Warning: skipping malformed alncl line: {line}\n")
                continue
            entries.append((parts[0].strip(), parts[1].strip(), ''))
    return entries

def extract_family_from_mirna(mirna):
    """Extract family number from a miRNA name like 'Ach-miR162a' -> '162'"""
    m = re.search(r'miR(\d+)', mirna)
    return m.group(1) if m else None

def parse_multi(filepath, start_num, known_families, prefix):
    """
    Process multi file, assign families based on majority voting.
    'known_families' contains numeric strings that are already taken (from -f FASTA).
    'prefix' is used only for newly assigned sequential numbers.
    For clusters that can reuse their intrinsic family number, the fixed prefix 'MIR' is used.
    """
    clusters = defaultdict(list)   # cluster_id -> list of (mirna, extra)
    cluster_order = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                sys.stderr.write(f"Warning: skipping malformed multi line: {line}\n")
                continue
            mirna = parts[0].strip()
            cluster_id = parts[1].strip()
            extra = parts[2].strip() if len(parts) >= 3 else ''
            if cluster_id not in clusters:
                cluster_order.append(cluster_id)
            clusters[cluster_id].append((mirna, extra))

    def cluster_sort_key(cid):
        m = re.search(r'\d+', cid)
        return int(m.group()) if m else 0
    cluster_order.sort(key=cluster_sort_key)

    # Collect statistics
    cluster_info = []
    for cid in cluster_order:
        members = clusters[cid]
        fam_counts = Counter()
        for mirna, _ in members:
            fam = extract_family_from_mirna(mirna)
            if fam:
                fam_counts[fam] += 1
        if fam_counts:
            max_count = max(fam_counts.values())
            best_fams = sorted([fam for fam, cnt in fam_counts.items() if cnt == max_count],
                               key=lambda x: int(x) if x.isdigit() else -1)
            best_fam = best_fams[-1]   # largest number on tie
        else:
            best_fam = None
        cluster_info.append({
            'cid': cid,
            'best_fam': best_fam,
            'size': len(members),
            'members': members
        })

    # Eligible clusters can use their intrinsic family if not in known_families
    eligible = [c for c in cluster_info if c['best_fam'] is not None
                and c['best_fam'] not in known_families]
    eligible.sort(key=lambda x: (-x['size'], cluster_sort_key(x['cid'])))

    assigned_fam = set()
    cluster_family_map = {}

    for c in eligible:
        fam = c['best_fam']
        if fam not in assigned_fam:
            # Always use MIR prefix for reused family numbers
            cluster_family_map[c['cid']] = f"MIR{fam}"
            assigned_fam.add(fam)

    # Remaining clusters get sequential numbers with the given prefix
    next_num = start_num
    used_numbers = set()
    for c in cluster_info:
        if c['cid'] not in cluster_family_map:
            fam_str = f"{prefix}{next_num}"
            cluster_family_map[c['cid']] = fam_str
            next_num += 1
        # record numeric parts (strip prefix)
        fam_str = cluster_family_map[c['cid']]
        m = re.search(r'MIRN?(\d+)', fam_str)
        if m:
            used_numbers.add(m.group(1))

    entries = []
    for c in cluster_info:
        fam = cluster_family_map[c['cid']]
        for mirna, extra in c['members']:
            entries.append((mirna, fam, extra))

    return entries, next_num, used_numbers

def parse_single(filepath, start_num, known_families, used_families, prefix):
    """
    Process single file. Each miRNA tries to reuse its intrinsic family if
    not in known_families and not already used. If possible, use prefix 'MIR';
    otherwise fall back to sequential numbering with the given prefix.
    """
    entries = []
    mirnas = []
    with open(filepath, 'r') as f:
        for line in f:
            mirna = line.strip()
            if not mirna:
                continue
            mirnas.append(mirna)

    used_single_fams = set(used_families)
    current_num = start_num

    for mirna in mirnas:
        fam = extract_family_from_mirna(mirna)
        if fam and fam not in known_families and fam not in used_single_fams:
            # Reuse family number with fixed MIR prefix
            family = f"MIR{fam}"
            used_single_fams.add(fam)
        else:
            # Assign a new number using the selected prefix
            family = f"{prefix}{current_num}"
            current_num += 1
        entries.append((mirna, family, ''))

    return entries

def main():
    parser = argparse.ArgumentParser(
        description='Create miRNA-to-family mapping with intelligent naming.')
    parser.add_argument('--alncl', default=None, help='Alncl file (two columns: mirna, family)')
    parser.add_argument('--single', default=None, help='Single file (one miRNA per line)')
    parser.add_argument('--multi', default=None, help='Multi (cluster) file (two or three columns)')
    parser.add_argument('--type', choices=['MIR', 'MIRN'], default='MIR',
                        help='Prefix for newly assigned families (default: MIR)')
    parser.add_argument('-s', '--start', type=int, default=None,
                        help='Starting number for sequential families (required when --type MIR)')
    parser.add_argument('-f', '--fasta', default=None,
                        help='FASTA file with known miRNA families')
    parser.add_argument('-o', '--output', required=True, help='Output mapping file')
    args = parser.parse_args()

    # At least one input file must be provided
    if not (args.alncl or args.single or args.multi):
        sys.exit("Error: at least one of --alncl, --single, --multi must be provided.")

    # Determine prefix and start number based on --type
    if args.type == 'MIR':
        prefix = 'MIR'
        if args.start is None:
            sys.exit("Error: -s/--start is required when --type MIR is used.")
        start_num = args.start
    else:  # MIRN
        prefix = 'MIRN'
        if args.start is not None:
            sys.exit("Error: -s/--start cannot be used with --type MIRN.")
        start_num = 1

    known_families = set()
    if args.fasta:
        known_families = parse_fasta_families(args.fasta)
        sys.stderr.write(f"Loaded {len(known_families)} known families from {args.fasta}\n")

    all_entries = []
    used_families = set()

    # 1. alncl (unchanged)
    if args.alncl:
        all_entries.extend(parse_alncl(args.alncl))

    # 2. multi
    next_num_for_single = start_num
    if args.multi:
        multi_entries, next_num_for_single, used_families = parse_multi(
            args.multi, start_num, known_families, prefix)
        all_entries.extend(multi_entries)

    # 3. single
    if args.single:
        single_entries = parse_single(args.single, next_num_for_single,
                                      known_families, used_families, prefix)
        all_entries.extend(single_entries)

    with open(args.output, 'w') as out:
        for mirna, family, extra in all_entries:
            out.write(f"{mirna}\t{family}\t{extra}\n")

    print(f"Output written to {args.output}", file=sys.stderr)

if __name__ == '__main__':
    main()
