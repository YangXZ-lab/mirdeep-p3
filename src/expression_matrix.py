#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
Build miRNA expression matrix from alignment files.

For each sample, parse the mature miRNA alignment file, filter alignments that
start within 2 bp of the 5' end and cover within 2 bp of the 3' end, calculate
both RPM (Reads Per Million) and raw read counts for each miRNA, and output
two matrices: one with RPM values and one with raw counts.

Inputs:
  -aln: one or more alignment files (comma-separated)
  -c: total read counts file (sample<TAB>total_reads)
  -species: species name (spaces replaced by underscores, used in output header)
  -prefix: comma-separated list of sample identifiers (must match -aln count)
  -o: output matrix file (RPM values)
  -oc: output matrix file (raw counts)
"""

import sys
import re
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional


def parse_total_reads(counts_file: str) -> Dict[str, int]:
    """Parse total read counts file."""
    total_reads = {}
    try:
        with open(counts_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    print(f"Warning: Line {line_num} has <2 columns, skipping.", file=sys.stderr)
                    continue
                sample, reads_str = parts[0], parts[1]
                try:
                    total_reads[sample] = int(reads_str)
                except ValueError:
                    print(f"Warning: Invalid total reads '{reads_str}' for sample '{sample}', skipping.", file=sys.stderr)
    except Exception as e:
        print(f"Error reading {counts_file}: {e}", file=sys.stderr)
        sys.exit(1)
    return total_reads


def extract_sample_name(aln_file: str) -> str:
    """Extract sample name from alignment file when -prefix not provided."""
    base = Path(aln_file).stem
    for suffix in ['.fa.mature', '.mature']:
        if base.endswith(suffix):
            base = base[:-len(suffix)]
            break
    if not base:
        base = Path(aln_file).name.replace('.aln', '')
    return base


def parse_alignment_file(aln_file: str, total_reads: int) -> Tuple[Dict[Tuple[str, str], float],
                                                                    Dict[Tuple[str, str], int]]:
    """
    Parse an alignment file and compute RPM and raw counts for each miRNA.

    Returns:
        rpm_dict: { (mir_name, mir_acc): total_rpm }
        count_dict: { (mir_name, mir_acc): total_raw_counts }
    """
    rpm_dict = {}
    count_dict = {}
    try:
        with open(aln_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                fields = line.split('\t')
                if len(fields) < 5:
                    print(f"Warning: Line {line_num} has <5 fields, skipping.", file=sys.stderr)
                    continue

                read_id = fields[0]          # e.g., read00072658_x4
                mir_info = fields[2]         # e.g., "23=N710-MIR390a=MIR_ID000000034"
                offset_str = fields[3]       # offset (integer)
                seq = fields[4]              # aligned sequence

                # Extract raw read count from read_id
                read_match = re.search(r'_x(\d+)', read_id)
                if not read_match:
                    continue
                raw_reads = int(read_match.group(1))

                # Extract miRNA length from mir_info (e.g., "23=")
                len_match = re.match(r'^(\d+)=', mir_info)
                if not len_match:
                    continue
                mir_len = int(len_match.group(1))

                # Remove length prefix to get "N710-MIR390a=MIR_ID000000034"
                mir_id_acc = re.sub(r'^\d+=', '', mir_info)
                # Split into name and accession
                if '=' in mir_id_acc:
                    mir_name, mir_acc = mir_id_acc.split('=', 1)
                else:
                    mir_name = mir_id_acc
                    mir_acc = ''

                # Parse offset
                try:
                    offset = int(offset_str)
                except ValueError:
                    continue

                # Filtering: offset > 2 or offset + len(seq) < mir_len - 2
                if offset > 2 or (offset + len(seq) < mir_len - 2):
                    continue

                key = (mir_name, mir_acc)
                rpm = raw_reads / total_reads * 1_000_000
                rpm_dict[key] = rpm_dict.get(key, 0.0) + rpm
                count_dict[key] = count_dict.get(key, 0) + raw_reads

    except Exception as e:
        print(f"Error processing {aln_file}: {e}", file=sys.stderr)
        sys.exit(1)

    return rpm_dict, count_dict


def write_matrix(output_file: str, sample_names: List[str],
                 all_keys: List[Tuple[str, str]],
                 sample_data: Dict[str, Dict[Tuple[str, str], float]],
                 value_type: str = 'rpm'):
    """
    Write a matrix (RPM or raw counts) to output file.

    Args:
        output_file: Path to output file.
        sample_names: Sorted list of sample names.
        all_keys: Sorted list of (mir_name, mir_acc) tuples.
        sample_data: Dict[sample][key] -> value (float for RPM, int for counts)
        value_type: 'rpm' or 'count' (for header comment)
    """
    try:
        with open(output_file, 'w') as fh:
            # Header: MIR, MIR_ID, then sample names
            header = ['MIR', 'MIR_ID'] + sample_names
            fh.write('\t'.join(header) + '\n')
            for mir_name, mir_acc in all_keys:
                row = [mir_name, mir_acc]
                for sample in sample_names:
                    val = sample_data[sample].get((mir_name, mir_acc), 0.0)
                    if value_type == 'rpm':
                        row.append(str(val))
                    else:  # raw counts (integer)
                        row.append(str(int(val)))
                fh.write('\t'.join(row) + '\n')
    except Exception as e:
        print(f"Error writing {output_file}: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Build miRNA expression matrices (RPM and raw counts)."
    )
    parser.add_argument('-aln', required=True,
                        help="Comma-separated list of alignment files.")
    parser.add_argument('-c', required=True,
                        help="Total read counts file (sample<TAB>total_reads).")
    parser.add_argument('-species', required=True,
                        help="Species name (spaces replaced by underscores).")
    parser.add_argument('-prefix', default=None,
                        help="Comma-separated list of sample identifiers (must match -aln count).")
    parser.add_argument('-o', required=True,
                        help="Output RPM matrix file.")
    parser.add_argument('-oc', required=True,
                        help="Output raw count matrix file.")
    args = parser.parse_args()

    # Process species
    species = args.species.replace(' ', '_')

    # Parse total read counts
    total_reads_dict = parse_total_reads(args.c)

    # Parse alignment files list
    aln_files = [f.strip() for f in args.aln.split(',') if f.strip()]
    if not aln_files:
        print("Error: No alignment files provided.", file=sys.stderr)
        sys.exit(1)

    # Determine sample names
    if args.prefix:
        sample_names = [p.strip() for p in args.prefix.split(',') if p.strip()]
        if len(sample_names) != len(aln_files):
            print(f"Error: Number of prefixes ({len(sample_names)}) does not match number of alignment files ({len(aln_files)}).", file=sys.stderr)
            sys.exit(1)
    else:
        sample_names = [extract_sample_name(f) for f in aln_files]

    # Validate sample names exist in counts file
    missing = [s for s in sample_names if s not in total_reads_dict]
    if missing:
        print(f"Error: Sample(s) not found in counts file: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)

    # Store data
    rpm_data = {sample: {} for sample in sample_names}
    count_data = {sample: {} for sample in sample_names}
    all_keys = set()

    for aln_file, sample in zip(aln_files, sample_names):
        total = total_reads_dict[sample]
        print(f"Processing {sample} from {aln_file} (total reads: {total:,})", file=sys.stderr)
        rpm_dict, count_dict = parse_alignment_file(aln_file, total)
        rpm_data[sample] = rpm_dict
        count_data[sample] = count_dict
        all_keys.update(rpm_dict.keys())

    if not all_keys:
        print("Error: No valid miRNA data found.", file=sys.stderr)
        sys.exit(1)

    # Sort keys and samples
    sorted_keys = sorted(all_keys, key=lambda x: (x[0], x[1]))
    sorted_samples = sorted(sample_names)

    # Write RPM matrix
    write_matrix(args.o, sorted_samples, sorted_keys, rpm_data, value_type='rpm')
    # Write raw count matrix
    write_matrix(args.oc, sorted_samples, sorted_keys, count_data, value_type='count')

    print(f"RPM matrix written to {args.o}", file=sys.stderr)
    print(f"Raw count matrix written to {args.oc}", file=sys.stderr)
    print(f"  Samples: {len(sorted_samples)}", file=sys.stderr)
    print(f"  miRNAs: {len(sorted_keys)}", file=sys.stderr)


if __name__ == "__main__":
    main()
