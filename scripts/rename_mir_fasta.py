#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
rename_mir_fasta.py

Read a FASTA file where each sequence header has the format:
    >MIR<old_number>-isoform<isoform_num>-<suffix>
and convert each header to:
    >MIR<counter>-isoform1-<suffix>
where <counter> starts from a user‑supplied value and increments by 1 for each sequence.
The sequence lines are copied unchanged.

Usage:
    python rename_mir_fasta.py -i input.fasta -s 11100 -o output.fasta
"""

import sys
import re
import argparse

def parse_header(header):
    """
    Parse a FASTA header line (without the leading '>').
    Expects format: MIR<number>-isoform<iso>-<suffix>
    Returns a tuple (old_number, isoform_num, suffix) if successful, else None.
    """
    header = header.strip()
    # Pattern: MIR followed by digits, then -isoform, then digits, then - and any suffix
    pattern = r'^MIR(\d+)-isoform(\d+)-(.+)$'
    match = re.match(pattern, header)
    if not match:
        return None
    old_number = int(match.group(1))
    isoform_num = int(match.group(2))
    suffix = match.group(3)
    return (old_number, isoform_num, suffix)

def generate_new_header(old_number, isoform_num, suffix, counter):
    """
    Generate the new header line (without the leading '>')
    using the current counter.
    Format: MIR<counter>-isoform1-<suffix>
    """
    return f"MIR{counter}-isoform1-{suffix}"

def main():
    parser = argparse.ArgumentParser(
        description="Rename MIR isoform FASTA headers by replacing the numeric part with a running counter."
    )
    parser.add_argument('-i', '--input', required=True,
                        help='Input FASTA file')
    parser.add_argument('-s', '--start', required=True, type=int,
                        help='Starting number for replacement (e.g., 11100)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output FASTA file')
    args = parser.parse_args()

    try:
        with open(args.input, 'r') as infile, open(args.output, 'w') as outfile:
            counter = args.start
            for line in infile:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    # Process header
                    header_content = line[1:]  # remove leading '>'
                    parsed = parse_header(header_content)
                    if parsed is None:
                        # Keep original header if it doesn't match the expected pattern
                        sys.stderr.write(f"Warning: header does not match expected pattern, keeping unchanged: {line}\n")
                        outfile.write(line + '\n')
                    else:
                        old_number, isoform_num, suffix = parsed
                        new_header = generate_new_header(old_number, isoform_num, suffix, counter)
                        outfile.write(f">{new_header}\n")
                        counter += 1
                else:
                    # Sequence line: copy unchanged
                    outfile.write(line + '\n')
    except IOError as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
