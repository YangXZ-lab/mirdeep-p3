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

Modified to match the behavior of the original Perl script exactly.
"""

import sys
import os
import re
import argparse

def parse_alignment_file(alignment_file):
    """
    Parse bowtie alignment file and extract read IDs that have alignments
    
    Args:
        alignment_file: Path to bowtie alignment output file
        
    Returns:
        Dictionary of read IDs that have alignments (value is always "T" to match Perl)
    """
    aligned_reads = {}
    
    try:
        with open(alignment_file, 'r') as f:
            for line in f:
                line = line.rstrip('\n\r')  # Match Perl's chomp and s/\r//
                if not line:
                    continue
                # Perl splits by tab and takes first element: my @a=split "\t"; $$in_hash{$a[0]}="T"
                parts = line.split('\t')
                if parts:
                    read_id = parts[0]  # Take entire first column, don't split by space
                    aligned_reads[read_id] = "T"  # Match Perl's value assignment
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
            reads_dict: Dictionary of {read_id: sequence} for filtered reads
            total_reads_count: Total number of reads (counting ALL reads, including ncRNA - matches Perl)
    """
    reads_dict = {}
    total_reads = 0
    current_id = None
    current_seq = []
    
    try:
        with open(reads_file, 'r') as f:
            for line in f:
                line = line.rstrip('\n\r')  # Match Perl's chomp and s/\r//
                if not line:
                    continue
                    
                if line.startswith('>'):
                    # Process previous read if exists
                    if current_id is not None:
                        sequence = ''.join(current_seq).upper()
                        seq_len = len(sequence)
                        
                        # Check if read is ncRNA - if yes, skip adding to reads_dict
                        # But note: total_reads was already incremented when we first saw this ID
                        if ncRNA_reads.get(current_id) != "T":
                            # Filter by length
                            if min_len <= seq_len <= max_len:
                                reads_dict[current_id] = sequence
                        
                    # Start new read
                    # Perl: chomp; s/\r//; s/^>//; $id=$_;
                    # So we remove '>' and keep everything else
                    current_id = line[1:]  # Remove '>' and keep entire line
                    current_seq = []
                    
                    # Extract count from ID (format: read00001_x123)
                    # Perl does this later in the main loop, but we need it here for total_reads
                    count_match = re.search(r'_x(\d+)$', current_id)
                    if count_match:
                        read_count = int(count_match.group(1))
                        total_reads += read_count  # Count ALL reads, including ncRNA (matches Perl)
                else:
                    current_seq.append(line)
            
            # Process the last read in file
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
    
    # Extract count from read ID
    count_match = re.search(r'_x(\d+)$', read_id)
    if count_match:
        count = int(count_match.group(1))
        rpm = (count / total_reads) * 1000000
        return rpm
    else:
        # If no count found, assume count = 1
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
                # Write specific reads
                for read_id in read_ids:
                    if read_id in reads_dict:
                        f.write(f">{read_id}\n{reads_dict[read_id]}\n")
            else:
                # Write all reads
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
    
    # Step 1: Parse alignment files to get ncRNA and known miRNA reads
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
    processed_reads = []  # All non-ncRNA reads for signature preparation
    filtered_reads = []   # Reads for precursor excision (high RPM or known miRNA)
    
    print("4. Filtering reads by RPM threshold...")
    high_rpm_count = 0
    known_miR_count = 0
    
    for read_id in reads_dict.keys():
        # Add to processed reads (all non-ncRNA reads)
        processed_reads.append(read_id)
        
        # Check if read maps to known miRNA
        # Note: Perl uses eq "T", so we check if the key exists with value "T"
        if known_miR_reads.get(read_id) == "T":
            filtered_reads.append(read_id)
            known_miR_count += 1
        else:
            # Check RPM threshold - match Perl's calculation exactly
            # Perl: $1/$total_reads*1000000 >= $threshold
            count_match = re.search(r'_x(\d+)$', read_id)
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
    
    # Write processed reads (all non-ncRNA reads)
    # Perl writes all reads from %reads_hash with ">$id\n$reads_hash{$id}\n"
    write_fasta(reads_dict, args.processed_file)
    print(f"   Processed reads written to: {args.processed_file}")
    
    # Write filtered reads (for precursor excision)
    write_fasta(reads_dict, args.filtered_file, filtered_reads)
    print(f"   Filtered reads written to: {args.filtered_file}")
    
    # Write total reads count - must match Perl exactly
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
    
    # Calculate approximate RPM distribution
    if filtered_reads:
        sample_rpms = []
        for i, read_id in enumerate(filtered_reads[:10]):  # Sample first 10
            count_match = re.search(r'_x(\d+)$', read_id)
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
