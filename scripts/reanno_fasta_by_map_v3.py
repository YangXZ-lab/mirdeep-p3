#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
rename_fasta_by_mapping.py

Rename FASTA sequence headers based on a two- or three-column mapping file.

New naming rule:
  - The target family (second column of the mapping) is used as the base name.
  - Within each family, sequences are assigned letter suffixes in the order
    they appear in the FASTA file: a, b, c, ..., z, aa, ab, ...
  - If --prefix is provided, the new name becomes 'prefix-familyletter'.
    Otherwise it is just 'familyletter'.

Input mapping file format (tab-separated, at least 2 columns):
    original_name   target_family   [extra]

Example:
    Amo-miR156o     MIR157
    =>  >MIR157a

    with --prefix "Ath" :
    =>  >Ath-MIR157a

Additional -r option writes an updated mapping file where original names are
replaced with the new names generated during renaming.

Usage:
    python rename_fasta_by_mapping.py -m mapping.txt -i input.fasta -o output.fasta [-r updated_map.txt] [--prefix "Ath"]
"""

import sys
import argparse
import re
from collections import defaultdict

def parse_mapping(mapping_file):
    """
    Read mapping file (tab‑separated, at least 2 columns).
    Returns:
        mapping: dict {original_name: target_family}
        lines: list of raw lines (list of fields) for later output
    """
    mapping = {}
    lines = []
    try:
        with open(mapping_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    sys.stderr.write(f"Warning: line {line_num} has <2 columns, skipping: {line}\n")
                    continue
                orig = parts[0].strip()
                target = parts[1].strip()
                mapping[orig] = target
                lines.append(parts)   # keep all columns
    except IOError as e:
        sys.stderr.write(f"Error reading mapping file: {e}\n")
        sys.exit(1)
    return mapping, lines

def num_to_letters(n):
    """Convert 0-based index to letter sequence: 0->'a', 1->'b', ..., 25->'z', 26->'aa', ..."""
    result = ''
    while True:
        n, rem = divmod(n, 26)
        result = chr(ord('a') + rem) + result
        if n == 0:
            break
        n -= 1
    return result

def process_fasta(input_fasta, mapping, output_handle, orig_to_new=None, prefix=None):
    """
    Read FASTA, rename headers according to mapping and family order,
    write to output_handle.
    If orig_to_new dict is provided, it will be filled with {orig_name: new_header}.
    """
    family_counter = defaultdict(int)   # family -> next letter index (0-based)
    
    current_header = None
    current_seq_lines = []
    seq_count = 0
    
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
                if current_header is not None:
                    seq_count += 1
                    new_header = generate_new_header(
                        current_header, mapping, family_counter, prefix
                    )
                    output_handle.write(f">{new_header}\n")
                    output_handle.write(''.join(current_seq_lines) + '\n')
                    if orig_to_new is not None:
                        orig_to_new[current_header] = new_header
                
                current_header = line[1:].strip()
                current_seq_lines = []
            else:
                if current_header is not None:
                    current_seq_lines.append(line)
        
        if current_header is not None:
            seq_count += 1
            new_header = generate_new_header(
                current_header, mapping, family_counter, prefix
            )
            output_handle.write(f">{new_header}\n")
            output_handle.write(''.join(current_seq_lines) + '\n')
            if orig_to_new is not None:
                orig_to_new[current_header] = new_header
    finally:
        if close_input:
            f.close()
    
    sys.stderr.write(f"Processed {seq_count} sequences.\n")

def generate_new_header(orig_name, mapping, counter_dict, prefix=None):
    """
    Generate new header based on mapping and family order.
    """
    if orig_name not in mapping:
        sys.stderr.write(f"Warning: '{orig_name}' not found in mapping. Keeping original name.\n")
        return orig_name
    
    target_family = mapping[orig_name]
    # Get the current index for this family, then assign the corresponding letter
    idx = counter_dict[target_family]
    letter = num_to_letters(idx)
    counter_dict[target_family] += 1
    
    if prefix:
        new_name = f"{prefix}-{target_family}{letter}"
    else:
        new_name = f"{target_family}{letter}"
    return new_name

def main():
    parser = argparse.ArgumentParser(
        description='Rename FASTA headers using a mapping file and family order.')
    parser.add_argument('-m', '--mapping', required=True,
                        help='Mapping file (original_name, target_family [, extra])')
    parser.add_argument('-i', '--input', help='Input FASTA file (default: stdin)')
    parser.add_argument('-o', '--output', help='Output FASTA file (default: stdout)')
    parser.add_argument('-r', '--rename-map', default=None,
                        help='Output updated mapping file with new names')
    parser.add_argument('--prefix', default=None,
                        help='Optional prefix for new names (e.g., "Ath")')
    args = parser.parse_args()
    
    # Read mapping (both dict and raw lines)
    mapping, mapping_lines = parse_mapping(args.mapping)
    sys.stderr.write(f"Loaded {len(mapping)} mapping entries.\n")
    
    # Prepare input / output streams
    if args.input:
        try:
            infile = open(args.input, 'r')
        except IOError as e:
            sys.stderr.write(f"Error opening input FASTA: {e}\n")
            sys.exit(1)
    else:
        infile = sys.stdin
    
    if args.output:
        try:
            outfile = open(args.output, 'w')
        except IOError as e:
            sys.stderr.write(f"Error opening output file: {e}\n")
            sys.exit(1)
    else:
        outfile = sys.stdout
    
    # Dictionary to collect old -> new name mappings
    orig_to_new = {} if args.rename_map else None
    
    # Process FASTA
    process_fasta(infile, mapping, outfile, orig_to_new, args.prefix)
    
    # Close files
    if infile is not sys.stdin:
        infile.close()
    if outfile is not sys.stdout:
        outfile.close()
    
    # Write updated mapping file if requested
    if args.rename_map:
        with open(args.rename_map, 'w') as f:
            for parts in mapping_lines:
                orig = parts[0].strip()
                # Replace original name with new name if available
                if orig in orig_to_new:
                    parts[0] = orig_to_new[orig]
                else:
                    sys.stderr.write(f"Note: original name '{orig}' not found in FASTA, keeping original.\n")
                f.write('\t'.join(parts) + '\n')
        sys.stderr.write(f"Updated mapping written to {args.rename_map}\n")

if __name__ == "__main__":
    main()
