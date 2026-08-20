#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
convert_bowtie_to_blast.py - Convert Bowtie alignment output to BLAST-like format

This script converts Bowtie output format to a BLAST-like tabular format.
It processes three input files:
  1. Bowtie alignment result file
  2. Short sequence file (reads) in FASTA format
  3. Chromosome/contig reference file in FASTA format

The output is in a BLAST-parsed format suitable for downstream processing.
"""

import sys
import os
import re
from collections import defaultdict
from typing import Dict, Tuple
import argparse

def parse_fasta_lengths(fasta_file: str) -> Dict[str, int]:
    """
    Parse FASTA file and extract sequence IDs and their lengths.
    
    Args:
        fasta_file: Path to FASTA file
        
    Returns:
        Dictionary mapping sequence ID to sequence length
    """
    seq_lengths = {}
    current_id = None
    current_seq = []
    
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_id is not None:
                        sequence = ''.join(current_seq)
                        seq_lengths[current_id] = len(sequence)
                    
                    # Start new sequence
                    # Extract ID (first word after '>')
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Save the last sequence
            if current_id is not None and current_seq:
                sequence = ''.join(current_seq)
                seq_lengths[current_id] = len(sequence)
        
        return seq_lengths
        
    except Exception as e:
        sys.stderr.write(f"Error: Failed to parse FASTA file {fasta_file}: {e}\n")
        sys.exit(1)

def parse_bowtie_line(line: str) -> Tuple[str, str, str, int, str]:
    """
    Parse a single line of Bowtie output.
    
    Bowtie format (tab-separated):
    read_id, strand, chromosome, position, sequence, quality, [optional fields]
    
    Args:
        line: A line from Bowtie output
        
    Returns:
        Tuple: (read_id, strand, chromosome, position, sequence)
    """
    try:
        fields = line.strip().split('\t')
        
        if len(fields) < 4:
            raise ValueError(f"Invalid Bowtie line: insufficient fields ({len(fields)} < 4)")
        
        read_id = fields[0].split()[0]  # Take first word only
        strand = fields[1]
        chromosome = fields[2]
        position = int(fields[3])
        sequence = fields[4] if len(fields) > 4 else ""
        
        return read_id, strand, chromosome, position, sequence
        
    except Exception as e:
        sys.stderr.write(f"Error: Failed to parse Bowtie line: {line}\n{str(e)}\n")
        raise

def calculate_alignment_coordinates(read_id: str, strand: str, position: int, 
                                  read_length: int, chr_length: int) -> Tuple[int, int, str]:
    """
    Calculate the alignment coordinates on the chromosome.
    
    Args:
        read_id: Read ID
        strand: Strand (+ or -)
        position: Alignment position from Bowtie
        read_length: Length of the read
        chr_length: Length of the chromosome
        
    Returns:
        Tuple: (alignment_start, alignment_end, strand_info) in 1-based coordinates
    """
    if strand == '+':
        # For plus strand: position is 0-based offset of first matching base
        alignment_start = position + 1  # Convert to 1-based
        alignment_end = position + read_length
        strand_info = "Plus/Plus"
        
    elif strand == '-':
        # For minus strand: position is 0-based offset of the 5' most base of the alignment
        # We use the original Bowtie position directly (common practice)
        alignment_start = position + 1  # Convert to 1-based (same as plus strand)
        alignment_end = position + read_length
        strand_info = "Plus/Minus"
        
        # Note: For negative strand, the read aligns to the reverse complement
        # but the coordinates are reported on the positive strand
        # This is the common practice in many alignment tools
    else:
        sys.stderr.write(f"Error: Invalid strand '{strand}' for read {read_id}\n")
        sys.exit(1)
    
    return alignment_start, alignment_end, strand_info

def validate_alignment_coordinates(start: int, end: int, chr_length: int, read_id: str, chromosome: str) -> bool:
    """
    Validate that alignment coordinates are within chromosome bounds.
    
    Args:
        start: Alignment start position (1-based)
        end: Alignment end position (1-based)
        chr_length: Chromosome length
        read_id: Read ID for error message
        chromosome: Chromosome ID for error message
        
    Returns:
        True if coordinates are valid, False otherwise
    """
    if start < 1:
        sys.stderr.write(f"Warning: Alignment start ({start}) < 1 for read {read_id} on {chromosome}\n")
        return False
    if end > chr_length:
        sys.stderr.write(f"Warning: Alignment end ({end}) > chromosome length ({chr_length}) for read {read_id} on {chromosome}\n")
        return False
    if start > end:
        sys.stderr.write(f"Warning: Alignment start ({start}) > end ({end}) for read {read_id} on {chromosome}\n")
        return False
    return True

def convert_bowtie_to_blast(bowtie_file: str, reads_lengths: Dict[str, int], 
                           chr_lengths: Dict[str, int]) -> str:
    """
    Convert Bowtie alignment file to BLAST-like format.
    
    Args:
        bowtie_file: Path to Bowtie output file
        reads_lengths: Dictionary of read IDs to their lengths
        chr_lengths: Dictionary of chromosome IDs to their lengths
        
    Returns:
        String containing the BLAST-formatted output
    """
    output_lines = []
    processed_count = 0
    skipped_count = 0
    warning_count = 0
    
    try:
        with open(bowtie_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                
                try:
                    # Parse Bowtie line
                    read_id, strand, chromosome, position, sequence = parse_bowtie_line(line)
                    
                    # Get read length
                    read_length = reads_lengths.get(read_id)
                    if read_length is None:
                        sys.stderr.write(f"Warning: Read ID '{read_id}' not found in reads file. Skipping.\n")
                        skipped_count += 1
                        continue
                    
                    # Get chromosome length
                    chr_length = chr_lengths.get(chromosome)
                    if chr_length is None:
                        sys.stderr.write(f"Warning: Chromosome '{chromosome}' not found in reference. Skipping.\n")
                        skipped_count += 1
                        continue
                    
                    # Calculate alignment coordinates
                    alignment_start, alignment_end, strand_info = calculate_alignment_coordinates(
                        read_id, strand, position, read_length, chr_length
                    )
                    
                    # Validate coordinates
                    if not validate_alignment_coordinates(alignment_start, alignment_end, 
                                                         chr_length, read_id, chromosome):
                        warning_count += 1
                    
                    # Format output in BLAST-like format
                    # Fields: query_id, query_length, query_range, subject_id, subject_length,
                    #         subject_range, e-value, identity, bitscore, strand_info
                    query_range = f"1..{read_length}"
                    subject_range = f"{alignment_start}..{alignment_end}"
                    
                    # Fixed values (as in original script)
                    e_value = "1e-04"
                    identity = "1.00"
                    bitscore = "42.1"
                    
                    # Create output line
                    output_line = (
                        f"{read_id}\t{read_length}\t{query_range}\t"
                        f"{chromosome}\t{chr_length}\t{subject_range}\t"
                        f"{e_value}\t{identity}\t{bitscore}\t{strand_info}"
                    )
                    
                    output_lines.append(output_line)
                    processed_count += 1
                    
                    # Progress feedback for large files
                    if line_num % 100000 == 0:
                        sys.stderr.write(f"  Processed {line_num:,} lines...\n")
                        
                except Exception as e:
                    sys.stderr.write(f"Warning: Skipping line {line_num}: {str(e)}\n")
                    skipped_count += 1
                    continue
        
        return "\n".join(output_lines), processed_count, skipped_count, warning_count
        
    except Exception as e:
        sys.stderr.write(f"Error: Failed to read Bowtie file {bowtie_file}: {e}\n")
        sys.exit(1)

def write_output(output_text: str, output_file: str = None):
    """
    Write output to file or stdout.
    
    Args:
        output_text: Text to write
        output_file: Output file path (None for stdout)
    """
    try:
        if output_file:
            with open(output_file, 'w') as f:
                f.write(output_text)
            print(f"Output written to: {output_file}")
        else:
            sys.stdout.write(output_text)
            
    except Exception as e:
        sys.stderr.write(f"Error: Failed to write output: {e}\n")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='Convert Bowtie alignment output to BLAST-like format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Input files:
  bowtie_result:  Bowtie alignment output file
  short_seq:      FASTA file of short reads (with sequences)
  chromosome:     FASTA file of chromosome/reference sequences
  
Output format (tab-separated):
  query_id, query_length, query_range, subject_id, subject_length,
  subject_range, e-value, identity, bitscore, strand_info
  
Important Note:
  For negative strand alignments, the coordinates are reported on the reference
  positive strand (common practice). The strand_info field indicates the actual
  strand orientation of the alignment.
  
Example Bowtie input line (positive strand):
  AtFlower100010_x2    +    MIR319c    508    AAGGAGATTCTTTCAGTCCAG
  
Example Bowtie input line (negative strand):
  AtFlower100011_x2    -    MIR319c    508    AAGGAGATTCTTTCAGTCCAG

Example output line (positive strand):
  AtFlower100010_x2    21    1..21    MIR319c    1000    509..529    1e-04    1.00    42.1    Plus/Plus

Example output line (negative strand):
  AtFlower100011_x2    21    1..21    MIR319c    1000    509..529    1e-04    1.00    42.1    Plus/Minus
  
Note: Negative strand alignments show coordinates on the positive strand reference.
      The actual alignment is to the reverse complement strand.
        """
    )
    
    parser.add_argument('bowtie_result',
                       help='Bowtie alignment output file')
    parser.add_argument('short_seq',
                       help='FASTA file of short reads')
    parser.add_argument('chromosome',
                       help='FASTA file of chromosome sequences')
    
    parser.add_argument('-o', '--output',
                       help='Output file (default: stdout)')
    
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Show detailed progress information')
    
    parser.add_argument('--skip-warnings', action='store_true',
                       help='Skip alignments with coordinate warnings')
    
    args = parser.parse_args()
    
    print("Starting conversion from Bowtie to BLAST format")
    print(f"Input files:")
    print(f"  Bowtie result: {args.bowtie_result}")
    print(f"  Short reads:   {args.short_seq}")
    print(f"  Chromosomes:   {args.chromosome}")
    if args.output:
        print(f"  Output file:   {args.output}")
    else:
        print(f"  Output:        stdout")
    print("-" * 60)
    
    # Step 1: Parse sequence lengths
    print("1. Parsing short reads file...")
    reads_lengths = parse_fasta_lengths(args.short_seq)
    print(f"   Found {len(reads_lengths):,} unique reads")
    
    print("2. Parsing chromosome file...")
    chr_lengths = parse_fasta_lengths(args.chromosome)
    print(f"   Found {len(chr_lengths):,} chromosomes/contigs")
    
    # Step 2: Convert Bowtie format
    print("3. Converting Bowtie alignments to BLAST format...")
    if args.verbose:
        print(f"   Reading {args.bowtie_result}...")
        print(f"   Note: Negative strand coordinates are reported on the reference positive strand")
    
    output_text, processed_count, skipped_count, warning_count = convert_bowtie_to_blast(
        args.bowtie_result, reads_lengths, chr_lengths
    )
    
    # Step 3: Write output
    print("4. Writing output...")
    write_output(output_text, args.output)
    
    # Step 4: Summary
    print("\n" + "=" * 60)
    print("CONVERSION COMPLETE - SUMMARY")
    print("=" * 60)
    print(f"Total alignments processed: {processed_count:,}")
    print(f"Alignments skipped:         {skipped_count:,}")
    print(f"Coordinate warnings:        {warning_count:,}")
    
    if processed_count > 0:
        success_rate = processed_count / (processed_count + skipped_count) * 100
        print(f"Success rate:               {success_rate:.1f}%")
    
    # Show sample of output
    if processed_count > 0 and args.verbose:
        print(f"\nSample output (first 3 lines):")
        lines = output_text.strip().split('\n')[:3]
        for i, line in enumerate(lines, 1):
            print(f"  {i}. {line}")
    
    print("=" * 60)

if __name__ == "__main__":
    main()
