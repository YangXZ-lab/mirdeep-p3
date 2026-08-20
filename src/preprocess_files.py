#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
Preprocessing module: Convert FastQ files to FastA format with count information
Main functions:
1. Read FastQ files
2. Merge duplicate sequences and count occurrences
3. Output FastA format with ID format: >read00001_x123
"""

import sys
import os
import gzip
from collections import defaultdict
from typing import Dict, List, Tuple
import argparse

def parse_fastq(file_path: str) -> Dict[str, int]:
    """
    Parse FastQ file and count sequence occurrences
    
    Args:
        file_path: Path to FastQ file
        
    Returns:
        Dictionary: sequence -> occurrence count
    """
    sequence_counts = defaultdict(int)
    
    # Check if file is gzipped
    if file_path.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'  # Text mode
    else:
        open_func = open
        mode = 'r'
    
    try:
        with open_func(file_path, mode) as f:
            line_num = 0
            sequence = None
            
            for line in f:
                line = line.strip()
                line_num += 1
                
                if line_num % 4 == 1:
                    continue
                elif line_num % 4 == 2:
                    sequence = line.upper()
                elif line_num % 4 == 0:
                    if sequence:
                        sequence_counts[sequence] += 1
                    sequence = None
    
    except Exception as e:
        print(f"Error: Cannot read file {file_path}: {e}", file=sys.stderr)
        sys.exit(1)
    
    return dict(sequence_counts)

def parse_fasta(file_path: str) -> Dict[str, int]:
    """
    Parse existing FastA file (if already in count format)
    
    Args:
        file_path: Path to FastA file
        
    Returns:
        Dictionary: sequence -> occurrence count
    """
    sequence_counts = {}
    
    try:
        with open(file_path, 'r') as f:
            current_seq = None
            
            for line in f:
                line = line.strip()
                
                if line.startswith('>'):
                    # Parse count information (format: >read00001_x123)
                    if '_x' in line:
                        try:
                            count = int(line.split('_x')[-1].split()[0])
                            seq_line = next(f).strip()
                            sequence_counts[seq_line.upper()] = count
                        except (ValueError, StopIteration):
                            print(f"Warning: Cannot parse count information: {line}", file=sys.stderr)
                elif line:
                    if current_seq is None:
                        current_seq = line.upper()
                        sequence_counts[current_seq] = 1
                    else:
                        current_seq += line.upper()
    
    except Exception as e:
        print(f"Error: Cannot read FastA file {file_path}: {e}", file=sys.stderr)
        sys.exit(1)
    
    return sequence_counts

def write_fasta_with_counts(sequence_counts: Dict[str, int], 
                           output_file: str, 
                           prefix: str = "read") -> Tuple[int, int]:
    """
    Write sequence count information to FastA file
    
    Args:
        sequence_counts: Dictionary of sequence->count
        output_file: Output file path
        prefix: Sequence ID prefix
        
    Returns:
        (total unique sequences, total reads)
    """
    total_reads = sum(sequence_counts.values())
    unique_sequences = len(sequence_counts)
    
    try:
        with open(output_file, 'w') as f:
            for i, (seq, count) in enumerate(sorted(sequence_counts.items()), 1):
                # Format: >read00001_x123
                seq_id = f">{prefix}{i:08d}_x{count}"
                f.write(f"{seq_id}\n{seq}\n")
    
    except Exception as e:
        print(f"Error: Cannot write output file {output_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    return unique_sequences, total_reads

def process_single_file(input_file: str, output_dir: str, 
                       input_type: str = 'fastq',
                       prefix: str = "read") -> Tuple[str, int, int]:
    """
    Process single input file
    
    Args:
        input_file: Input file path
        output_dir: Output directory
        input_type: Input file type ('fastq' or 'fasta')
        prefix: Sequence ID prefix
        
    Returns:
        (output file path, unique sequences count, total reads count)
    """
    # Extract filename (without extension)
    filename = os.path.splitext(os.path.basename(input_file))[0]
    if filename.endswith('.fq') or filename.endswith('.fastq'):
        filename = os.path.splitext(filename)[0]
    
    # Create output filename
    output_file = os.path.join(output_dir, f"{filename}.processed.fa")
    
    # Parse input file
    if input_type.lower() == 'fastq':
        sequence_counts = parse_fastq(input_file)
    else:  # fasta
        sequence_counts = parse_fasta(input_file)
    
    # Write output file
    unique_seqs, total_reads = write_fasta_with_counts(
        sequence_counts, output_file, prefix
    )
    
    return output_file, unique_seqs, total_reads

def process_multiple_files(file_list: List[str], output_dir: str,
                          input_type: str = 'fastq',
                          prefix: str = "read") -> Dict[str, dict]:
    """
    Process multiple input files
    
    Args:
        file_list: List of input file paths
        output_dir: Output directory
        input_type: Input file type
        prefix: Sequence ID prefix
        
    Returns:
        Dictionary: filename -> processing information
    """
    results = {}
    
    for input_file in file_list:
        if not os.path.exists(input_file):
            print(f"Warning: File does not exist, skipping: {input_file}", file=sys.stderr)
            continue
        
        output_file, unique_seqs, total_reads = process_single_file(
            input_file, output_dir, input_type, prefix
        )
        
        results[input_file] = {
            'output_file': output_file,
            'unique_sequences': unique_seqs,
            'total_reads': total_reads
        }
        
        print(f"Processing completed: {input_file}")
        print(f"  Output file: {output_file}")
        print(f"  Unique sequences: {unique_seqs:,}")
        print(f"  Total reads: {total_reads:,}")
        print()
    
    return results

def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Preprocess sRNA-seq data: Convert FastQ/FastA to FastA format with count information',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input file(s) (comma-separated) or text file containing file list'
    )
    
    parser.add_argument(
        '-o', '--output',
        default='.',
        help='Output directory (default: current directory)'
    )
    
    parser.add_argument(
        '-t', '--input-type',
        choices=['fastq', 'fasta'],
        default='fastq',
        help='Input file type: fastq or fasta (default: fastq)'
    )
    
    parser.add_argument(
        '-b', '--batch',
        action='store_true',
        help='Input is a batch file with one file path per line'
    )
    
    parser.add_argument(
        '--prefix',
        default='read',
        help='Sequence ID prefix (default: read)'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Prepare file list
    file_list = []
    
    if args.batch:
        try:
            with open(args.input, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        file_list.append(line)
        except Exception as e:
            print(f"Error: Cannot read batch file {args.input}: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        # Parse comma-separated string
        file_list = [f.strip() for f in args.input.split(',') if f.strip()]
    
    if not file_list:
        print("Error: No valid input files found", file=sys.stderr)
        sys.exit(1)
    
    # Process files
    print(f"Starting to process {len(file_list)} file(s)...")
    print(f"Input type: {args.input_type}")
    print(f"Output directory: {args.output}")
    print("-" * 50)
    
    results = process_multiple_files(
        file_list, args.output, args.input_type, args.prefix
    )
    
    # Output statistics
    print("\n" + "=" * 50)
    print("Processing completed! Summary statistics:")
    print("=" * 50)
    
    total_unique = 0
    total_reads = 0
    
    for input_file, info in results.items():
        total_unique += info['unique_sequences']
        total_reads += info['total_reads']
        print(f"{os.path.basename(input_file)}:")
        print(f"  Unique sequences: {info['unique_sequences']:,}")
        print(f"  Total reads: {info['total_reads']:,}")
    
    print(f"\nTotal:")
    print(f"  Total unique sequences: {total_unique:,}")
    print(f"  Total reads: {total_reads:,}")

if __name__ == "__main__":
    main()
