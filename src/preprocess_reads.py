#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.


"""
preprocess_reads.py - Filter sRNA-seq reads for miRNA prediction

This script preprocesses and filters original reads file for miRDeep-P pipeline.
It filters reads by:
1. Removing reads that map to ncRNA (rRNA, tRNA, snRNA)
2. Considering reads that map to known miRNA sequences
3. Applying RPM threshold to retain high-abundance reads

"""

import sys
import os
import re
import argparse

# Precompile regex once (avoids re-compilation in every loop iteration)
COUNT_RE = re.compile(r'_x(\d+)$')


def parse_alignment_file(alignment_file):
    """
    Parse bowtie alignment file and extract read IDs that have alignments

    Args:
        alignment_file: Path to bowtie alignment output file

    Returns:
        Dictionary of read IDs that have alignments (value is always "T")
    """
    aligned_reads = {}

    try:
        with open(alignment_file, 'r') as f:
            for line in f:
                line = line.rstrip('\n\r')
                if not line:
                    continue
                parts = line.split('\t')
                if parts:
                    aligned_reads[parts[0]] = "T"
        return aligned_reads
    except Exception as e:
        sys.stderr.write(f"Error: Failed to read alignment file {alignment_file}: {e}\n")
        sys.exit(1)


def parse_reads_file(reads_file, ncRNA_reads, min_len=19, max_len=24):
    """
    Parse reads file and filter by length and ncRNA mapping

    Args:
        reads_file: Path to reads file in fasta format
        ncRNA_reads: Dictionary of read IDs that map to ncRNA
        min_len: Minimum read length (default: 19)
        max_len: Maximum read length (default: 24)

    Returns:
        tuple: (reads_dict, total_reads_count)
    """
    reads_dict = {}
    total_reads = 0
    current_id = None
    current_seq = []

    try:
        with open(reads_file, 'r') as f:
            for line in f:
                line = line.rstrip('\n\r')
                if not line:
                    continue

                if line.startswith('>'):
                    if current_id is not None:
                        sequence = ''.join(current_seq).upper()
                        seq_len = len(sequence)

                        if ncRNA_reads.get(current_id) != "T":
                            if min_len <= seq_len <= max_len:
                                reads_dict[current_id] = sequence

                    current_id = line[1:]
                    current_seq = []

                    count_match = COUNT_RE.search(current_id)
                    if count_match:
                        total_reads += int(count_match.group(1))
                else:
                    current_seq.append(line)

            if current_id is not None and current_seq:
                sequence = ''.join(current_seq).upper()
                seq_len = len(sequence)

                if ncRNA_reads.get(current_id) != "T" and min_len <= seq_len <= max_len:
                    reads_dict[current_id] = sequence

        return reads_dict, total_reads

    except Exception as e:
        sys.stderr.write(f"Error: Failed to read file {reads_file}: {e}\n")
        sys.exit(1)


def calculate_rpm(read_id, total_reads):
    """
    Calculate RPM (Reads Per Million) for a read

    Args:
        read_id: Read ID string containing count (format: read00001_x123)
        total_reads: Total number of reads in the dataset

    Returns:
        RPM value as float
    """
    if total_reads == 0:
        return 0.0

    count_match = COUNT_RE.search(read_id)
    if count_match:
        count = int(count_match.group(1))
        rpm = (count / total_reads) * 1000000
        return rpm
    else:
        return (1 / total_reads) * 1000000


def write_fasta(reads_dict, output_file, read_ids=None):
    """
    Write reads to fasta file

    Args:
        reads_dict: Dictionary of {read_id: sequence}
        output_file: Path to output fasta file
        read_ids: Optional list of specific read IDs to write
    """
    try:
        with open(output_file, 'w') as f:
            if read_ids is not None:
                for read_id in read_ids:
                    if read_id in reads_dict:
                        f.write(f">{read_id}\n{reads_dict[read_id]}\n")
            else:
                for read_id, sequence in reads_dict.items():
                    f.write(f">{read_id}\n{sequence}\n")
    except Exception as e:
        sys.stderr.write(f"Error: Failed to write output file {output_file}: {e}\n")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Preprocess and filter sRNA-seq reads for miRNA prediction',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Input files:
  reads_file:     FastA format unique reads file with IDs like '>read00001_x123'
  rfam_aln:       Bowtie alignment output of reads against ncRNA sequences
  miR_aln:        Bowtie alignment output of reads against known miRNA sequences
  threshold:      RPM threshold for retaining reads for precursor excision

Output files:
  processed_file: Reads for signature preparation (all non-ncRNA reads)
  filtered_file:  Reads for precursor excision (high-RPM or known miRNA reads)
  total_reads:    Total reads count file
        """
    )

    parser.add_argument('reads_file',
                       help='FastA format reads file with unique sequences')
    parser.add_argument('rfam_aln',
                       help='Bowtie alignment file for ncRNA (rRNA/tRNA/snRNA)')
    parser.add_argument('miR_aln',
                       help='Bowtie alignment file for known miRNA sequences')
    parser.add_argument('threshold', type=float,
                       help='RPM threshold for reads to be retrieved for precursor excision')
    parser.add_argument('processed_file',
                       help='Output path for processed reads (for signature preparation)')
    parser.add_argument('filtered_file',
                       help='Output path for filtered reads (for precursor excision)')
    parser.add_argument('total_reads_file',
                       help='Output path for total reads count')

    parser.add_argument('--min-len', type=int, default=19,
                       help='Minimum read length (default: 19)')
    parser.add_argument('--max-len', type=int, default=24,
                       help='Maximum read length (default: 24)')

    args = parser.parse_args()

    print(f"Starting preprocessing of reads file: {args.reads_file}")
    print(f"Parameters:")
    print(f"  ncRNA alignment file: {args.rfam_aln}")
    print(f"  miRNA alignment file: {args.miR_aln}")
    print(f"  RPM threshold: {args.threshold}")
    print(f"  Read length range: {args.min_len}-{args.max_len} nt")
    print("-" * 60)

    # Step 1: Parse alignment files
    print("1. Parsing ncRNA alignment file...")
    ncRNA_reads = parse_alignment_file(args.rfam_aln)
    print(f"   Found {len(ncRNA_reads)} reads mapping to ncRNA")

    print("2. Parsing known miRNA alignment file...")
    known_miR_reads = parse_alignment_file(args.miR_aln)
    print(f"   Found {len(known_miR_reads)} reads mapping to known miRNA")

    # Step 2: Parse reads file and filter by length and ncRNA
    print("3. Parsing reads file and filtering by length...")
    reads_dict, total_reads = parse_reads_file(
        args.reads_file, ncRNA_reads, args.min_len, args.max_len
    )
    print(f"   Total reads (counting ALL reads, including ncRNA): {total_reads:,}")
    print(f"   Unique reads after filtering (non-ncRNA, length {args.min_len}-{args.max_len}): {len(reads_dict):,}")

    # Step 3: Prepare lists for output files
    processed_reads = []
    filtered_reads = []

    print("4. Filtering reads by RPM threshold...")
    high_rpm_count = 0
    known_miR_count = 0

    for read_id in reads_dict.keys():
        processed_reads.append(read_id)

        if known_miR_reads.get(read_id) == "T":
            filtered_reads.append(read_id)
            known_miR_count += 1
        else:
            count_match = COUNT_RE.search(read_id)
            if count_match:
                count = int(count_match.group(1))
                rpm = (count / total_reads) * 1000000
                if rpm >= args.threshold:
                    filtered_reads.append(read_id)
                    high_rpm_count += 1

    print(f"   Reads mapping to known miRNA: {known_miR_count}")
    print(f"   Reads meeting RPM threshold: {high_rpm_count}")
    print(f"   Total reads for precursor excision: {len(filtered_reads):,}")

    # Step 4: Write output files
    print("5. Writing output files...")

    write_fasta(reads_dict, args.processed_file)
    print(f"   Processed reads written to: {args.processed_file}")

    write_fasta(reads_dict, args.filtered_file, filtered_reads)
    print(f"   Filtered reads written to: {args.filtered_file}")

    try:
        with open(args.total_reads_file, 'w') as f:
            f.write(str(total_reads))
        print(f"   Total reads count written to: {args.total_reads_file}")
    except Exception as e:
        sys.stderr.write(f"Error: Failed to write total reads file: {e}\n")
        sys.exit(1)

    # Step 5: Calculate and report statistics
    filtered_percentage = (len(filtered_reads) / len(reads_dict)) * 100 if reads_dict else 0

    print("\n" + "=" * 60)
    print("PROCESSING COMPLETE - SUMMARY")
    print("=" * 60)
    print(f"Total reads in dataset (including ncRNA): {total_reads:,}")
    print(f"Unique non-ncRNA reads (after length filtering): {len(reads_dict):,}")
    print(f"Reads for signature preparation: {len(processed_reads):,}")
    print(f"Reads for precursor excision: {len(filtered_reads):,}")
    print(f"Filtered reads percentage: {filtered_percentage:.1f}%")

    if filtered_reads:
        sample_rpms = []
        for i, read_id in enumerate(filtered_reads[:10]):
            count_match = COUNT_RE.search(read_id)
            if count_match:
                count = int(count_match.group(1))
                rpm = (count / total_reads) * 1000000
                sample_rpms.append(rpm)

        if sample_rpms:
            avg_rpm = sum(sample_rpms) / len(sample_rpms)
            print(f"Average RPM of filtered reads (sample): {avg_rpm:.1f}")

    print("=" * 60)


if __name__ == "__main__":
    main()
