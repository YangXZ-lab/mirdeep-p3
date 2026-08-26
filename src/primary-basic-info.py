#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
Generate miRNA primary information (including miRNA* and extended precursor regions).

Inputs:
  -f filtered-all-nr-anno (tab-separated, 11+ columns)
  -struc_20nt stemloop_20nt.struc (RNAfold output for 20nt extension)
  -struc stemloop.struc (RNAfold output for full precursor)
  -chr_length chromosome length file (two columns: chrom<TAB>length)
  -species species name (spaces will be replaced by underscores)
  -prefix_miRNA prefix for miRNA IDs (spaces replaced by underscores)
  -m mapping file (optional, for MIRN family assignment)
  -o output file path

Output: tab-separated file with 24 columns (ranges split into start/end).
"""

import sys
import re
import argparse
from typing import Dict, Tuple, Optional
from collections import defaultdict


def parse_range(range_str: str) -> Tuple[int, int]:
    """
    Parse a range string like "start..end" into (start, end) integers.
    """
    parts = range_str.split('..')
    if len(parts) != 2:
        raise ValueError(f"Invalid range format: {range_str}")
    return int(parts[0]), int(parts[1])


def parse_struc_file(file_path: str, remove_strand: bool = True) -> Dict[str, Dict[str, str]]:
    """
    Parse RNAfold output file (three-line records: >id, sequence, structure).

    Args:
        file_path: Path to the .struc file.
        remove_strand: If True, strip trailing '(+)/'(-)' from the ID.

    Returns:
        Dictionary: id -> {'seq': sequence, 'struc': structure}
    """
    data = {}
    try:
        with open(file_path, 'r') as f:
            lines = [line.rstrip('\n') for line in f if line.strip()]
        i = 0
        while i < len(lines):
            if lines[i].startswith('>'):
                header = lines[i][1:].strip()
                if remove_strand:
                    header = re.sub(r'\([+-]\)$', '', header).strip()
                if i + 2 >= len(lines):
                    break
                seq = lines[i+1].strip().upper().replace('U', 'T')
                struc_line = lines[i+2].strip()
                struc = re.sub(r'\s*\(\s*[\d\.\-]+\s*\)$', '', struc_line).rstrip()
                data[header] = {'seq': seq, 'struc': struc}
                i += 3
            else:
                i += 1
    except Exception as e:
        print(f"Error reading {file_path}: {e}", file=sys.stderr)
        sys.exit(1)
    return data


def read_chr_lengths(file_path: str) -> Dict[str, int]:
    """
    Read chromosome length file (two columns: chrom<TAB>length).
    """
    lengths = {}
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    continue
                chrom = parts[0]
                try:
                    length = int(parts[1])
                    lengths[chrom] = length
                except ValueError:
                    print(f"Warning: Invalid length for {chrom}: {parts[1]}", file=sys.stderr)
    except Exception as e:
        print(f"Error reading {file_path}: {e}", file=sys.stderr)
        sys.exit(1)
    return lengths


def strip_strand(name: str) -> str:
    """Remove trailing strand annotation like (+), (-), etc."""
    return re.sub(r'\([+-]\)$', '', name).strip()


def load_name_to_family_map(map_file: str) -> Dict[str, str]:
    """
    Load a mapping file (two or three columns: name_with_strand, family, [weight]).
    The name is normalized by stripping the strand suffix.
    Returns dict: normalized_name -> family.
    """
    mapping = {}
    try:
        with open(map_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    print(f"Warning: Line {line_num} has <2 columns, skipping.", file=sys.stderr)
                    continue
                name_raw = parts[0].strip()
                family = parts[1].strip()
                name = strip_strand(name_raw)
                mapping[name] = family
    except Exception as e:
        print(f"Error reading mapping file {map_file}: {e}", file=sys.stderr)
        sys.exit(1)
    return mapping


def read_filtered_anno(file_path: str) -> Tuple[Dict[str, Dict[str, str]], Dict[str, str]]:
    """
    Read filtered-all-nr-anno file.

    Expected columns (0-indexed):
        0: chrom
        1: strand
        2: read_id
        3: name
        4: mature_range
        5: pre_range
        6: mature_seq
        7: pre_seq
        8: conservation (e.g., "non_conserved" or "conserved=MIR408")
        9: confidence (optional)
        10: source (optional)

    Returns:
        (conserve_dict, non_conserve_dict)
        conserve_dict: {family_number: {name: full_line}}
        non_conserve_dict: {name: full_line}
    """
    conserve = {}
    non_conserve = {}
    with open(file_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.rstrip('\n')
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) < 9:
                print(f"Warning: Line {line_num} has <9 columns, skipping.", file=sys.stderr)
                continue
            name = fields[3]
            conservation = fields[8]
            if conservation == "non_conserved":
                non_conserve[name] = line
            else:
                match = re.search(r'conserved=(MIR[N]?\d+)', conservation)
                if match:
                    fam = match.group(1)
                    num_match = re.search(r'(\d+)', fam)
                    fam_num = int(num_match.group(1)) if num_match else 0
                    if fam_num not in conserve:
                        conserve[fam_num] = {}
                    conserve[fam_num][name] = line
                else:
                    print(f"Warning: Cannot parse conservation '{conservation}' at line {line_num}", file=sys.stderr)
    return conserve, non_conserve


def generate_tag(member: int) -> str:
    """
    Generate letter suffix for miRNA member (a,b,c,...,z,aa,ab,...).
    member starts from 1.
    """
    member -= 1
    if member < 26:
        return chr(ord('a') + member)
    else:
        first = member // 26 - 1
        second = member % 26
        return chr(ord('a') + first) + chr(ord('a') + second)


def extract_family_number(family_str: str) -> Optional[int]:
    """Extract the numeric part from a family string like MIRN5 -> 5."""
    m = re.search(r'MIRN?(\d+)', family_str)
    return int(m.group(1)) if m else None


def process_conserved_records(conserve, stem_20, stem_full, chr_lengths,
                              species, prefix_miRNA, output_handle, accession_counter):
    """Process all conserved miRNA families."""
    for fam_num in sorted(conserve.keys()):
        any_name = next(iter(conserve[fam_num].keys()))
        any_line = conserve[fam_num][any_name]
        any_fields = any_line.split('\t')
        conservation = any_fields[8]
        match = re.search(r'conserved=(MIR[N]?\d+)', conservation)
        if not match:
            continue
        family_name = match.group(1)
        member = 1
        for name in sorted(conserve[fam_num].keys()):
            tag = generate_tag(member)
            full_line = conserve[fam_num][name]
            fields = full_line.split('\t')
            while len(fields) < 11:
                fields.append('')
            chrom = fields[0]
            strand = fields[1]
            read_id = fields[2]
            mature_range_str = fields[4]
            pre_range_str = fields[5]
            mature_seq = fields[6].upper()
            pre_seq = fields[7].upper().replace('U', 'T')
            confidence = fields[9] if len(fields) > 9 else '.'
            source = fields[10] if len(fields) > 10 else '.'

            try:
                mature_start, mature_end = parse_range(mature_range_str)
                pre_start, pre_end = parse_range(pre_range_str)
            except ValueError as e:
                print(f"Warning: {e} for name {name}, skipping.", file=sys.stderr)
                member += 1
                continue

            ext_start = max(1, pre_start - 20)
            ext_end = min(chr_lengths.get(chrom, pre_end + 20), pre_end + 20)

            stem_20_rec = stem_20.get(name, {'seq': '', 'struc': ''})
            stem_full_rec = stem_full.get(name, {'seq': '', 'struc': ''})

            ind = pre_seq.find(mature_seq)
            star_seq = ''
            star_beg = 0
            star_end = 0
            if ind >= 0 and abs(ind) <= 2:
                star_seq = pre_seq[-21:] if len(pre_seq) >= 21 else pre_seq
                star_beg = max(1, pre_end - 20)
                star_end = pre_end
            elif abs(ind - (len(pre_seq) - len(mature_seq) - 1)) <= 2:
                star_seq = pre_seq[:21] if len(pre_seq) >= 21 else pre_seq
                star_beg = pre_start + 1
                star_end = pre_start + 21
            else:
                print(f"Warning: Cannot determine star location for {name} (ind={ind}, pre_len={len(pre_seq)}, mature_len={len(mature_seq)}), skipping.", file=sys.stderr)
                member += 1
                continue

            star_seq = star_seq.upper()
            star_beg = max(1, star_beg)
            star_end = min(chr_lengths.get(chrom, star_end), star_end)

            mature_family = family_name.replace('MIR', 'miR')
            miR_id = f"{prefix_miRNA}-{family_name}{tag}"
            mature_mir_id = f"{prefix_miRNA}-{mature_family}{tag}"
            star_id = f"{prefix_miRNA}-{mature_family}{tag}*"

            out_fields = [
                miR_id, f"MIR_ID{accession_counter:09d}", family_name, species,
                chrom, str(ext_start), str(ext_end), strand,
                stem_20_rec['seq'], stem_20_rec['struc'],
                pre_seq, stem_full_rec['struc'],
                str(pre_start), str(pre_end),
                mature_mir_id, str(mature_start), str(mature_end), mature_seq,
                star_id, star_seq, str(star_beg), str(star_end),
                confidence, source
            ]
            output_handle.write('\t'.join(out_fields) + '\n')
            accession_counter += 1
            member += 1
    return accession_counter


def process_non_conserved_records(non_conserve, stem_20, stem_full, chr_lengths,
                                  species, prefix_miRNA, output_handle,
                                  accession_counter, name_to_family=None):
    """
    Process non-conserved records.
    If name_to_family is provided, use it to assign MIRN families:
    - For records whose name maps to a MIRN family, group by that family and
      assign letter suffixes in sorted name order.
    - Remaining records get sequential MIRN numbers, avoiding conflicts.
    """
    matched_family_to_names = defaultdict(list)   # family -> list of (name, line)
    unmatched_records = []

    for name, line in non_conserve.items():
        mapped_family = None
        if name_to_family:
            norm_name = strip_strand(name)
            mapped_family = name_to_family.get(norm_name)
        if mapped_family and mapped_family.startswith('MIRN'):
            matched_family_to_names[mapped_family].append((name, line))
        else:
            unmatched_records.append((name, line))

    # Collect used family numbers to avoid conflicts
    used_family_numbers = set()
    for fam in matched_family_to_names.keys():
        num = extract_family_number(fam)
        if num is not None:
            used_family_numbers.add(num)

    # Determine next available MIRN number for unmatched
    next_unmatched_num = 1
    while next_unmatched_num in used_family_numbers:
        next_unmatched_num += 1

    # Helper to process a single record
    def process_record(name, line, family_name, tag):
        nonlocal accession_counter
        fields = line.split('\t')
        while len(fields) < 11:
            fields.append('')
        chrom = fields[0]
        strand = fields[1]
        read_id = fields[2]
        mature_range_str = fields[4]
        pre_range_str = fields[5]
        mature_seq = fields[6].upper()
        pre_seq = fields[7].upper().replace('U', 'T')
        confidence = fields[9] if len(fields) > 9 else '.'
        source = fields[10] if len(fields) > 10 else '.'

        try:
            mature_start, mature_end = parse_range(mature_range_str)
            pre_start, pre_end = parse_range(pre_range_str)
        except ValueError as e:
            print(f"Warning: {e} for non-conserved {name}, skipping.", file=sys.stderr)
            return

        ext_start = max(1, pre_start - 20)
        ext_end = min(chr_lengths.get(chrom, pre_end + 20), pre_end + 20)

        stem_20_rec = stem_20.get(name, {'seq': '', 'struc': ''})
        stem_full_rec = stem_full.get(name, {'seq': '', 'struc': ''})

        ind = pre_seq.find(mature_seq)
        star_seq = ''
        star_beg = 0
        star_end = 0
        if ind != -1 and abs(ind) <= 2:
            star_seq = pre_seq[-21:] if len(pre_seq) >= 21 else pre_seq
            star_beg = max(1, pre_end - 20)
            star_end = pre_end
        elif ind != -1 and abs(ind - (len(pre_seq) - len(mature_seq) - 1)) <= 2:
            star_seq = pre_seq[:21] if len(pre_seq) >= 21 else pre_seq
            star_beg = pre_start + 1
            star_end = pre_start + 21
        else:
            print(f"Warning: Cannot determine star location for non-conserved {name}, skipping.", file=sys.stderr)
            return

        star_seq = star_seq.upper()
        star_beg = max(1, star_beg)
        star_end = min(chr_lengths.get(chrom, star_end), star_end)

        mature_family = family_name.replace('MIR', 'miR')
        miR_id = f"{prefix_miRNA}-{family_name}{tag}"
        mature_mir_id = f"{prefix_miRNA}-{mature_family}{tag}"
        star_id = f"{prefix_miRNA}-{mature_family}{tag}*"

        out_fields = [
            miR_id, f"MIR_ID{accession_counter:09d}", family_name, species,
            chrom, str(ext_start), str(ext_end), strand,
            stem_20_rec['seq'], stem_20_rec['struc'],
            pre_seq, stem_full_rec['struc'],
            str(pre_start), str(pre_end),
            mature_mir_id, str(mature_start), str(mature_end), mature_seq,
            star_id, star_seq, str(star_beg), str(star_end),
            confidence, source
        ]
        output_handle.write('\t'.join(out_fields) + '\n')
        accession_counter += 1

    # Process matched families (sorted by numeric family id)
    for fam in sorted(matched_family_to_names.keys(), key=extract_family_number):
        fam_records = sorted(matched_family_to_names[fam], key=lambda x: x[0])
        for idx, (name, line) in enumerate(fam_records):
            tag = generate_tag(idx + 1)   # a, b, c, ...
            process_record(name, line, fam, tag)

    # Process unmatched records (sorted by name)
    unmatched_records.sort(key=lambda x: x[0])
    for name, line in unmatched_records:
        fam = f"MIRN{next_unmatched_num}"
        next_unmatched_num += 1
        while next_unmatched_num in used_family_numbers:
            next_unmatched_num += 1
        process_record(name, line, fam, '')   # no suffix for single-member families

    return accession_counter


def main():
    parser = argparse.ArgumentParser(
        description="Generate miRNA primary information (including miRNA* and extended precursor regions).")
    parser.add_argument('-f', '--filtered-anno', required=True,
                        help="filtered-all-nr-anno file (tab-separated, 11+ columns)")
    parser.add_argument('-struc_20nt', required=True,
                        help="stemloop_20nt.struc file (RNAfold output)")
    parser.add_argument('-struc', required=True,
                        help="stemloop.struc file (RNAfold output)")
    parser.add_argument('-chr_length', required=True,
                        help="chromosome length file (two columns: chrom<TAB>length)")
    parser.add_argument('-species', required=True,
                        help="Species name (spaces will be replaced by underscores)")
    parser.add_argument('-prefix_miRNA', required=True,
                        help="Prefix for miRNA IDs (spaces replaced by underscores)")
    parser.add_argument('-m', '--map', default=None,
                        help="Mapping file for MIRN family assignment (optional)")
    parser.add_argument('-o', '--output', required=True,
                        help="Output file path")

    args = parser.parse_args()

    species = args.species.replace(' ', '_')
    prefix_miRNA = args.prefix_miRNA.replace(' ', '_')

    # Read optional mapping file
    name_to_family = None
    if args.map:
        name_to_family = load_name_to_family_map(args.map)

    # Read input files
    chr_lengths = read_chr_lengths(args.chr_length)
    stem_20 = parse_struc_file(args.struc_20nt, remove_strand=True)
    stem_full = parse_struc_file(args.struc, remove_strand=True)
    conserve, non_conserve = read_filtered_anno(args.filtered_anno)

    # Open output file
    try:
        out_fh = open(args.output, 'w')
    except Exception as e:
        print(f"Error: Cannot open output file {args.output}: {e}", file=sys.stderr)
        sys.exit(1)

    accession_counter = 1
    with out_fh:
        accession_counter = process_conserved_records(
            conserve, stem_20, stem_full, chr_lengths, species, prefix_miRNA,
            out_fh, accession_counter
        )
        accession_counter = process_non_conserved_records(
            non_conserve, stem_20, stem_full, chr_lengths, species, prefix_miRNA,
            out_fh, accession_counter, name_to_family
        )

    print(f"Processing completed. Output written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()