#!/usr/bin/env python3
"""
Combine cluster, singleton and alignment results into a three-column
miRNA-to-family mapping file.

Third column (ratio) is preserved for entries from the --multi file;
for other sources it is left empty.

New: -f FASTA file provides known families. Multi clusters are renamed
using the most frequent family (if tie, pick the largest numeric value)
if it is not already known and not conflicting.
Singletons are renamed similarly, with conflict resolution.

Usage:
    python merge_family.py --alncl alncl.txt --single single.txt --multi multi.txt \
                           -s 2001 -f known.fasta -o output.txt
"""

import argparse
import sys
import re
from collections import defaultdict, Counter

def parse_fasta_families(fasta_file):
    """Extract family identifiers from FASTA headers (>MIR2111... -> '2111')"""
    families = set()
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                m = re.match(r'>MIR?(\d+)', line)
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

def parse_multi(filepath, start_num, known_families):
    """
    Process multi file, assign families based on majority voting.
    For each cluster, the family with the most occurrences is chosen;
    if tie, the largest numeric family wins.
    The chosen family is used only if it is not in known_families and not
    already assigned to a larger or equal cluster (priority by cluster size, then CL order).
    Otherwise fall back to MIRN numbering.

    Returns:
        entries: list of (mirna, family, extra)
        next_num: int, next available MIRN number for singletons
        used_families: set of family numbers (strings) already used
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
            # Among ties, pick the largest numeric family
            best_fams = sorted([fam for fam, cnt in fam_counts.items() if cnt == max_count],
                               key=lambda x: int(x) if x.isdigit() else -1)
            best_fam = best_fams[-1]   # largest number
        else:
            best_fam = None
        cluster_info.append({
            'cid': cid,
            'best_fam': best_fam,
            'size': len(members),
            'members': members
        })

    # Eligibility: must have best_fam, not in known_families, and not outbid.
    eligible = [c for c in cluster_info if c['best_fam'] is not None
                and c['best_fam'] not in known_families]
    # Sort by cluster size descending, then by CL order (natural)
    eligible.sort(key=lambda x: (-x['size'], cluster_sort_key(x['cid'])))

    assigned_fam = set()
    cluster_family_map = {}

    for c in eligible:
        fam = c['best_fam']
        if fam not in assigned_fam:
            cluster_family_map[c['cid']] = f"MIR{fam}"
            assigned_fam.add(fam)

    # Remaining clusters receive MIRN numbers
    next_mirn = start_num
    used_mirn_numbers = set()
    for c in cluster_info:
        if c['cid'] not in cluster_family_map:
            fam_str = f"MIR{next_mirn}"
            cluster_family_map[c['cid']] = fam_str
            next_mirn += 1
        # Collect used numeric parts
        fam_str = cluster_family_map[c['cid']]
        m = re.search(r'MIR?(\d+)', fam_str)
        if m:
            used_mirn_numbers.add(m.group(1))

    entries = []
    for c in cluster_info:
        fam = cluster_family_map[c['cid']]
        for mirna, extra in c['members']:
            entries.append((mirna, fam, extra))

    return entries, next_mirn, used_mirn_numbers

def parse_single(filepath, start_num, known_families, used_families):
    """
    Process single file.
    Each miRNA tries to use its intrinsic family (from name) if not in
    known_families and not already used. Otherwise fall back to MIRN numbers.
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
    current_mirn = start_num

    for mirna in mirnas:
        fam = extract_family_from_mirna(mirna)
        if fam and fam not in known_families and fam not in used_single_fams:
            family = f"MIR{fam}"
            used_single_fams.add(fam)
        else:
            family = f"MIR{current_mirn}"
            current_mirn += 1
        entries.append((mirna, family, ''))

    return entries

def main():
    parser = argparse.ArgumentParser(
        description='Create miRNA-to-family mapping with intelligent naming.')
    parser.add_argument('--alncl', default=None, help='Alncl file (two columns: mirna, family)')
    parser.add_argument('--single', default=None, help='Single file (one miRNA per line)')
    parser.add_argument('--multi', default=None, help='Multi (cluster) file (two or three columns)')
    parser.add_argument('-s', '--start', type=int, required=True,
                        help='Starting MIR number for default numbering')
    parser.add_argument('-f', '--fasta', default=None,
                        help='FASTA file with known miRNA families')
    parser.add_argument('-o', '--output', required=True, help='Output mapping file')
    args = parser.parse_args()

    known_families = set()
    if args.fasta:
        known_families = parse_fasta_families(args.fasta)
        sys.stderr.write(f"Loaded {len(known_families)} known families from {args.fasta}\n")

    all_entries = []
    used_families = set()

    # 1. alncl
    if args.alncl:
        all_entries.extend(parse_alncl(args.alncl))

    # 2. multi
    next_num_for_single = args.start
    if args.multi:
        multi_entries, next_num_for_single, used_families = parse_multi(
            args.multi, args.start, known_families)
        all_entries.extend(multi_entries)

    # 3. single
    if args.single:
        single_entries = parse_single(args.single, next_num_for_single,
                                      known_families, used_families)
        all_entries.extend(single_entries)

    with open(args.output, 'w') as out:
        for mirna, family, extra in all_entries:
            out.write(f"{mirna}\t{family}\t{extra}\n")

    print(f"Output written to {args.output}", file=sys.stderr)

if __name__ == '__main__':
    main()
