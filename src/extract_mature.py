#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
Convert prediction file (8 columns) to mature miRNA BED file and extract sequences.

"""

import sys
import re
import argparse
import subprocess
import shutil
from typing import Optional, Tuple


def parse_range(range_str: str) -> Optional[Tuple[int, int]]:
    """
    Parse a range string "start..end" into (start, end) integers.

    Args:
        range_str: String like "148456427..148456448"

    Returns:
        Tuple (start, end) as integers, or None if format is invalid.
    """
    match = re.match(r'^(\d+)\.\.(\d+)$', range_str)
    if match:
        return int(match.group(1)), int(match.group(2))
    return None


def determine_signature(mature_beg: int, mature_end: int,
                        pre_beg: int, pre_end: int,
                        strand: str) -> str:
    """
    Determine the signature of mature miRNA relative to pre-miRNA.

    Args:
        mature_beg: Start coordinate of mature miRNA (original).
        mature_end: End coordinate of mature miRNA (original).
        pre_beg: Start coordinate of pre-miRNA (original).
        pre_end: End coordinate of pre-miRNA (original).
        strand: Strand direction ('+' or '-').

    Returns:
        '5' if mature aligns with 5' end of pre-miRNA,
        '3' if aligns with 3' end,
        'A' for ambiguous (internal) location.
    """
    if (mature_beg == pre_beg and strand == '+') or \
       (mature_end == pre_end and strand == '-'):
        return '5'
    elif (mature_beg == pre_beg and strand == '-') or \
         (mature_end == pre_end and strand == '+'):
        return '3'
    else:
        return 'A'


def process_line(fields: list) -> Optional[str]:
    """
    Process a single line from the input prediction file.

    Expected fields (0-indexed):
        0: chromosome
        1: strand
        2: read ID (with RPM)
        3: unique name
        4: mature range (start..end)
        5: pre range (start..end)
        6: mature sequence (unused)
        7: pre sequence (unused)

    Args:
        fields: List of 8 tab-separated fields.

    Returns:
        Formatted output line (6 columns, tab-separated) or None if line is invalid.
    """
    if len(fields) < 8:
        print(f"Warning: Skipping line with {len(fields)} columns (expected 8): {fields}",
              file=sys.stderr)
        return None

    chrom = fields[0]
    strand = fields[1]
    name = fields[3]
    mature_range_str = fields[4]
    pre_range_str = fields[5]

    # Parse mature range
    mature_coords = parse_range(mature_range_str)
    if mature_coords is None:
        print(f"Warning: Cannot parse mature range '{mature_range_str}', skipping line.",
              file=sys.stderr)
        return None
    mature_beg, mature_end = mature_coords

    # Parse pre range (needed for signature)
    pre_coords = parse_range(pre_range_str)
    if pre_coords is None:
        print(f"Warning: Cannot parse pre range '{pre_range_str}', skipping line.",
              file=sys.stderr)
        return None
    pre_beg, pre_end = pre_coords

    # Determine signature using original coordinates (no adjustments)
    sign = determine_signature(mature_beg, mature_end, pre_beg, pre_end, strand)

    # Output BED line: chrom, mature_start, mature_end, name, sign, strand
    out_line = f"{chrom}\t{mature_beg}\t{mature_end}\t{name}\t{sign}\t{strand}"
    return out_line


def convert_predictions_to_mature_bed(input_file: str, output_bed: str) -> Tuple[int, int]:
    """
    Convert the prediction file to a mature miRNA BED file.

    Args:
        input_file: Path to input file (8-column tab-separated).
        output_bed: Path to output BED file (6-column).

    Returns:
        Tuple (processed_lines, skipped_lines)
    """
    processed = 0
    skipped = 0

    try:
        with open(input_file, 'r') as infh, open(output_bed, 'w') as outfh:
            for line_num, line in enumerate(infh, 1):
                line = line.rstrip('\n')
                if not line:
                    continue

                fields = line.split('\t')
                out_line = process_line(fields)

                if out_line is None:
                    skipped += 1
                    continue

                outfh.write(out_line + '\n')
                processed += 1
    except Exception as e:
        print(f"Error: Failed to process file {input_file}: {e}", file=sys.stderr)
        sys.exit(1)

    return processed, skipped


def run_bedtools_getfasta(bed_file: str, genome_fasta: str, output_fasta: str) -> None:
    """
    Run bedtools getfasta to extract sequences from the genome.

    Args:
        bed_file: Input BED file (6 columns).
        genome_fasta: Genome FASTA file.
        output_fasta: Output FASTA file.

    Raises:
        SystemExit: If bedtools is not found or command fails.
    """
    # Check if bedtools is available
    bedtools_path = shutil.which('bedtools')
    if bedtools_path is None:
        print("Error: 'bedtools' not found in PATH. Please install bedtools or add it to PATH.",
              file=sys.stderr)
        sys.exit(1)

    cmd = [
        bedtools_path, 'getfasta',
        '-s',           # force strand (reverse complement for negative strand)
        '-nameOnly',    # use BED name column as FASTA header
        '-bed', bed_file,
        '-fi', genome_fasta,
        '-fo', output_fasta
    ]

    print(f"Running: {' '.join(cmd)}", file=sys.stderr)
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # bedtools prints progress to stderr; we can optionally show it
        if result.stderr:
            print(result.stderr, file=sys.stderr)
    except subprocess.CalledProcessError as e:
        print(f"Error: bedtools getfasta failed with exit code {e.returncode}", file=sys.stderr)
        print(f"stderr: {e.stderr}", file=sys.stderr)
        sys.exit(1)


def main():
    """Command line interface."""
    parser = argparse.ArgumentParser(
        description='Convert prediction file (8 columns) to mature miRNA BED file '
                    'and extract sequences using bedtools getfasta.',
        epilog='Coordinates are used as-is from the input (no adjustments). '
               'Output BED: chrom<TAB>mature_start<TAB>mature_end<TAB>name<TAB>signature<TAB>strand',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input file path (8-column tab-separated prediction file).'
    )

    parser.add_argument(
        '-o', '--output-bed',
        required=True,
        help='Output BED file path (6 columns).'
    )

    parser.add_argument(
        '-fo', '--output-fasta',
        required=True,
        help='Output FASTA file path (extracted mature miRNA sequences).'
    )

    parser.add_argument(
        '-genome', '--genome-fasta',
        required=True,
        help='Genome FASTA file (used by bedtools getfasta).'
    )

    args = parser.parse_args()

    # Step 1: Convert to BED
    print(f"Processing input file: {args.input}", file=sys.stderr)
    processed, skipped = convert_predictions_to_mature_bed(args.input, args.output_bed)

    print(f"Conversion completed.", file=sys.stderr)
    print(f"  Processed lines: {processed}", file=sys.stderr)
    if skipped:
        print(f"  Skipped lines: {skipped}", file=sys.stderr)
    print(f"  BED file written: {args.output_bed}", file=sys.stderr)

    if processed == 0:
        print("Error: No valid records found. Exiting without running bedtools.", file=sys.stderr)
        sys.exit(1)

    # Step 2: Extract sequences with bedtools
    run_bedtools_getfasta(args.output_bed, args.genome_fasta, args.output_fasta)
    print(f"FASTA file written: {args.output_fasta}", file=sys.stderr)
    print("All done.", file=sys.stderr)


if __name__ == "__main__":
    main()
