#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
rename_fasta_by_mapping.py

Rename FASTA sequence headers based on a mapping file and generate
sequential 'rename' numbers per target family.

Input mapping file format (tab/space separated, at least 3 columns):
    original_name   any   target_family

The original_name can be any string; the script will extract the
trailing numeric suffix (e.g., '-7', '-2') to preserve in the new name.

Output FASTA header format:
    >{target_family}-rename{seq_num}-{suffix}

If no suffix found, output: >{target_family}-rename{seq_num}

Usage:
    python rename_fasta_by_mapping.py -m mapping.txt -i input.fasta -o output.fasta
    cat input.fasta | python rename_fasta_by_mapping.py -m mapping.txt > output.fasta
"""

import sys
import argparse
import re
from collections import defaultdict

def parse_mapping(mapping_file):
    """
    Read mapping file (tab/space separated).
    Returns dict: {original_name: target_family}
    """
    mapping = {}
    try:
        with open(mapping_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) < 3:
                    sys.stderr.write(f"Warning: line {line_num} has <3 columns, skipping: {line}\n")
                    continue
                orig = parts[0]
                target = parts[2]   # third column
                mapping[orig] = target
    except IOError as e:
        sys.stderr.write(f"Error reading mapping file: {e}\n")
        sys.exit(1)
    return mapping

def extract_suffix(name):
    """
    Extract the trailing numeric suffix from a name.
    Looks for a hyphen followed by digits at the end of the string.
    Returns suffix as string (e.g., '7', '2').
    If no match, returns empty string.
    """
    # Match pattern: -(\d+)$
    match = re.search(r'-(\d+)$', name)
    if match:
        return match.group(1)
    else:
        # No suffix found, not an error, just return empty
        return ""

def process_fasta(input_fasta, mapping, output_handle):
    """
    Read FASTA, rename headers, write to output_handle.
    """
    family_counter = defaultdict(int)
    
    # State for reading FASTA
    current_header = None
    current_seq_lines = []
    seq_count = 0
    
    # Handle input source
    if hasattr(input_fasta, 'read'):
        f = input_fasta
        close_input = False
    else:
        f = open(input_fasta, 'r')
        close_input = True
    
    try:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                # Write previous sequence if any
                if current_header is not None:
                    seq_count += 1
                    new_header = generate_new_header(current_header, mapping, family_counter)
                    output_handle.write(f">{new_header}\n")
                    output_handle.write(''.join(current_seq_lines) + '\n')
                
                # Start new sequence
                current_header = line[1:].strip()  # remove '>'
                current_seq_lines = []
            else:
                if current_header is not None:
                    current_seq_lines.append(line)
                # else ignore lines before first '>'
        
        # Write last sequence
        if current_header is not None:
            seq_count += 1
            new_header = generate_new_header(current_header, mapping, family_counter)
            output_handle.write(f">{new_header}\n")
            output_handle.write(''.join(current_seq_lines) + '\n')
    finally:
        if close_input:
            f.close()
    
    sys.stderr.write(f"Processed {seq_count} sequences.\n")

def generate_new_header(orig_name, mapping, counter_dict):
    """
    Generate new header based on mapping and counter.
    """
    if orig_name not in mapping:
        sys.stderr.write(f"Warning: '{orig_name}' not found in mapping. Keeping original name.\n")
        return orig_name
    
    target_family = mapping[orig_name]
    suffix = extract_suffix(orig_name)
    
    # Increment counter for this family
    counter_dict[target_family] += 1
    seq_num = counter_dict[target_family]
    
    if suffix:
        return f"{target_family}-rename{seq_num}-{suffix}"
    else:
        return f"{target_family}-rename{seq_num}"

def main():
    parser = argparse.ArgumentParser(description='Rename FASTA headers using mapping file and sequential rename numbers.')
    parser.add_argument('-m', '--mapping', required=True, help='Mapping file (original_name, any, target_family)')
    parser.add_argument('-i', '--input', help='Input FASTA file (default: stdin)')
    parser.add_argument('-o', '--output', help='Output FASTA file (default: stdout)')
    args = parser.parse_args()
    
    # Read mapping
    mapping = parse_mapping(args.mapping)
    sys.stderr.write(f"Loaded {len(mapping)} mapping entries.\n")
    
    # Open input FASTA
    if args.input:
        try:
            infile = open(args.input, 'r')
        except IOError as e:
            sys.stderr.write(f"Error opening input FASTA: {e}\n")
            sys.exit(1)
    else:
        infile = sys.stdin
    
    # Open output
    if args.output:
        try:
            outfile = open(args.output, 'w')
        except IOError as e:
            sys.stderr.write(f"Error opening output file: {e}\n")
            sys.exit(1)
    else:
        outfile = sys.stdout
    
    # Process
    process_fasta(infile, mapping, outfile)
    
    # Cleanup
    if infile is not sys.stdin:
        infile.close()
    if outfile is not sys.stdout:
        outfile.close()

if __name__ == "__main__":
    main()
