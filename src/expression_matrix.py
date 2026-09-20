#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
Build miRNA expression matrix from alignment files.

For each sample, parse the mature miRNA alignment file, keep the alignments
that start within ``--max-5p-trunc`` bp of the 5' end and end within
``--max-3p-trunc`` bp of the 3' end, then compute both RPM (Reads Per Million)
and raw read counts for each miRNA.

When ``-fa`` is given, the reference FASTA used for the bowtie alignment
(``<prefix>-basic-info-mature-1nt.fa``) is parsed as well, so that every
reference miRNA receives a row in the output matrices even when no read in any
sample passes the coverage filter. Such rows are zero-filled, and the row order
follows the reference FASTA, which is the same order as
``<prefix>-basic-info``. Without ``-fa`` the row set is derived from the
alignments alone, reproducing the behaviour up to 3.1.6c.

Inputs:
  -aln: one or more alignment files (comma-separated)
  -c: total read counts file (sample<TAB>total_reads)
  -species: species name (spaces replaced by underscores, used in output header)
  -prefix: comma-separated list of sample identifiers (must match -aln count)
  -fa: reference FASTA used for the alignment (enables zero-filling)
  -o: output matrix file (RPM values)
  -oc: output matrix file (raw counts)
"""

import sys
import re
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional


# 3.1.6c used 5'<=2 and 3'<=2. The 3' tolerance is relaxed to 3 bp here because
# plant sRNA-seq libraries are dominated by 19-21 nt reads and a read that is
# only 1 bp shorter than the previous cut-off silently drops the whole miRNA
# out of the matrix.
DEFAULT_MAX_5P_TRUNC = 2
DEFAULT_MAX_3P_TRUNC = 3


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


def parse_reference_fasta(fa_file: str) -> Tuple[List[Tuple[str, str]], Dict[Tuple[str, str], int]]:
    """
    Parse the reference FASTA used for the alignment.

    Headers are written by ``extract_mature_flank.py`` as
    ``>{ref_len}={mir_name}={mir_acc}``. The returned key order follows the
    file, which in turn follows the row order of ``*-basic-info``.

    Returns:
        keys: ordered list of (mir_name, mir_acc)
        ref_len: { (mir_name, mir_acc): reference length }
    """
    keys: List[Tuple[str, str]] = []
    ref_len: Dict[Tuple[str, str], int] = {}
    seen = set()
    state = {'name': None, 'acc': '', 'declared': 0, 'seq_len': 0}

    def _flush():
        name = state['name']
        if name is None:
            return
        key = (name, state['acc'])
        ref_len[key] = state['seq_len'] or state['declared']
        if key not in seen:
            seen.add(key)
            keys.append(key)

    with open(fa_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                _flush()
                parts = line[1:].split('=', 2)
                if len(parts) >= 2:
                    try:
                        state['declared'] = int(parts[0])
                    except ValueError:
                        state['declared'] = 0
                    state['name'] = parts[1]
                    state['acc'] = parts[2] if len(parts) > 2 else ''
                else:
                    state['declared'] = 0
                    state['name'] = line[1:]
                    state['acc'] = ''
                state['seq_len'] = 0
            else:
                state['seq_len'] += len(line)
    _flush()
    return keys, ref_len


def parse_alignment_file(aln_file: str, total_reads: int,
                         max_5p_trunc: int = DEFAULT_MAX_5P_TRUNC,
                         max_3p_trunc: int = DEFAULT_MAX_3P_TRUNC,
                         ref_len: Optional[Dict[Tuple[str, str], int]] = None,
                         stats: Optional[Dict[str, int]] = None) -> Tuple[Dict[Tuple[str, str], float],
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

                # Extract reference length from mir_info (e.g., "23=")
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

                # Coverage filter. gap_5p is the number of reference bases in
                # front of the read, gap_3p the number left uncovered behind
                # it. A negative gap means the read runs past that end.
                gap_5p = offset
                gap_3p = mir_len - (offset + len(seq))
                if gap_5p > max_5p_trunc or gap_3p > max_3p_trunc:
                    if stats is not None:
                        stats['filtered'] = stats.get('filtered', 0) + 1
                        if gap_5p > max_5p_trunc:
                            stats['fail_5p'] = stats.get('fail_5p', 0) + 1
                        if gap_3p > max_3p_trunc:
                            stats['fail_3p'] = stats.get('fail_3p', 0) + 1
                    continue

                key = (mir_name, mir_acc)
                rpm = raw_reads / total_reads * 1_000_000
                rpm_dict[key] = rpm_dict.get(key, 0.0) + rpm
                count_dict[key] = count_dict.get(key, 0) + raw_reads
                if stats is not None:
                    stats['kept'] = stats.get('kept', 0) + 1

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
    parser.add_argument('-fa', '--reference-fasta', dest='fa', default=None,
                        help="Reference FASTA used for the alignment, i.e. "
                             "'<prefix>-basic-info-mature-1nt.fa'. When given, every "
                             "reference miRNA is written to the matrix (zero-filled if no "
                             "read passes the filter) in reference order.")
    parser.add_argument('--max-5p-trunc', type=int, default=DEFAULT_MAX_5P_TRUNC,
                        help="Maximum allowed 5' truncation in bp (default: %(default)s).")
    parser.add_argument('--max-3p-trunc', type=int, default=DEFAULT_MAX_3P_TRUNC,
                        help="Maximum allowed 3' truncation in bp (default: %(default)s; "
                             "use 2 to reproduce the 3.1.6c behaviour, 4 for the most "
                             "permissive setting tested).")
    parser.add_argument('-o', required=True,
                        help="Output RPM matrix file.")
    parser.add_argument('-oc', required=True,
                        help="Output raw count matrix file.")
    args = parser.parse_args()

    if args.max_5p_trunc < 0 or args.max_3p_trunc < 0:
        print("Error: truncation tolerances must be non-negative.", file=sys.stderr)
        sys.exit(1)

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

    # Optional reference FASTA for zero-filling
    ref_keys = None
    ref_len = None
    if args.fa:
        if not Path(args.fa).is_file():
            print(f"Error: reference FASTA not found: {args.fa}", file=sys.stderr)
            sys.exit(1)
        ref_keys, ref_len = parse_reference_fasta(args.fa)
        if not ref_keys:
            print(f"Error: no sequence found in reference FASTA: {args.fa}", file=sys.stderr)
            sys.exit(1)
        print(f"Reference miRNAs in {Path(args.fa).name}: {len(ref_keys)}", file=sys.stderr)

    # Store data
    rpm_data = {sample: {} for sample in sample_names}
    count_data = {sample: {} for sample in sample_names}
    all_keys = set()
    stats: Dict[str, int] = {}

    for aln_file, sample in zip(aln_files, sample_names):
        total = total_reads_dict[sample]
        print(f"Processing {sample} from {aln_file} (total reads: {total:,})", file=sys.stderr)
        rpm_dict, count_dict = parse_alignment_file(
            aln_file, total,
            max_5p_trunc=args.max_5p_trunc,
            max_3p_trunc=args.max_3p_trunc,
            ref_len=ref_len,
            stats=stats,
        )
        rpm_data[sample] = rpm_dict
        count_data[sample] = count_dict
        all_keys.update(rpm_dict.keys())

    if not all_keys:
        print("Error: No valid miRNA data found.", file=sys.stderr)
        sys.exit(1)

    kept = stats.get('kept', 0)
    filtered = stats.get('filtered', 0)
    total_lines = kept + filtered
    print(f"Coverage filter: 5'<={args.max_5p_trunc} bp, 3'<={args.max_3p_trunc} bp", file=sys.stderr)
    if total_lines:
        print(f"  alignments kept: {kept:,} / {total_lines:,} ({100.0 * kept / total_lines:.2f}%)", file=sys.stderr)
        print(f"  dropped: {filtered:,} (5' over: {stats.get('fail_5p', 0):,}, "
              f"3' over: {stats.get('fail_3p', 0):,}, categories may overlap)", file=sys.stderr)

    # Determine the row set. With -fa the reference defines it, so miRNAs that
    # only failed the coverage filter no longer disappear from the matrices.
    if ref_keys is not None:
        ordered_keys = list(ref_keys)
        extra = all_keys - set(ref_keys)
        if extra:
            print(f"Warning: {len(extra)} aligned miRNA(s) absent from the reference FASTA; appending them.",
                  file=sys.stderr)
            ordered_keys.extend(sorted(extra))
        zero_filled = len(set(ref_keys) - all_keys)
        print(f"  rows zero-filled: {zero_filled} (of {len(ref_keys)} reference miRNAs)", file=sys.stderr)
        sorted_keys = ordered_keys
    else:
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
