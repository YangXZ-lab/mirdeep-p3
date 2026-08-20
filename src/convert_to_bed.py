#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
Script to convert miRNA prediction file to BED format.
Function: Process miRNA prediction results and generate BED format files for visualization.
"""

import re
import argparse
import sys
from pathlib import Path

def process_line(line, line_num):
    """Process a single line of data"""
    line = line.rstrip('\n')
    fields = line.split('\t')
    
    if len(fields) < 8:
        print(f"Warning: Line {line_num} has insufficient fields, skipping", file=sys.stderr)
        return None
    
    # Parse mature and pre coordinates
    mature_match = re.match(r'^(\d+)\.\.(\d+)$', fields[4])
    pre_match = re.match(r'^(\d+)\.\.(\d+)$', fields[5])
    
    if not mature_match or not pre_match:
        print(f"Warning: Line {line_num} has invalid coordinate format, skipping", file=sys.stderr)
        return None
    
    mature_beg = int(mature_match.group(1))
    pre_beg = int(pre_match.group(1))
    
    # Calculate lengths and end positions
    mature_len = len(fields[6])
    pre_len = len(fields[7])
    
    mature_end = mature_beg + mature_len
    pre_end = pre_beg + pre_len
    
    # Determine sign type
    strand = fields[1]
    if (mature_beg == pre_beg and strand == "+") or (mature_end == pre_end and strand == "-"):
        sign = "5"
    elif (mature_beg == pre_beg and strand == "-") or (mature_end == pre_end and strand == "+"):
        sign = "3"
    else:
        sign = "A"
    
    # Return BED format line
    return f"{fields[0]}\t{pre_beg}\t{pre_end}\t{fields[3]}\t{sign}\t{strand}"

def main():
    parser = argparse.ArgumentParser(
        description='Convert miRNA prediction file to BED format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage examples:
  python convert_to_bed.py -i input_file -o output_file
  python convert_to_bed.py input_file  # Output to input_file.bed
        """
    )
    
    parser.add_argument('input_file', nargs='?', help='Input file path (filter_P_prediction)')
    parser.add_argument('-i', '--input', help='Input file path (filter_P_prediction)')
    parser.add_argument('-o', '--output', help='Output BED file path (default: input_file.bed)')
    
    args = parser.parse_args()
    
    # Determine input file
    input_file = args.input if args.input else args.input_file
    if not input_file:
        parser.print_help()
        print("\nError: Input file must be specified", file=sys.stderr)
        sys.exit(1)
    
    # Check if input file exists
    input_path = Path(input_file)
    if not input_path.exists():
        print(f"Error: Input file '{input_file}' does not exist", file=sys.stderr)
        sys.exit(1)
    
    # Determine output file
    if args.output:
        output_file = args.output
    else:
        output_file = str(input_path) + '.bed'
    
    # Process files
    processed_count = 0
    skipped_count = 0
    
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            line_num = 0
            for line in infile:
                line_num += 1
                
                # Skip empty lines
                if not line.strip():
                    continue
                
                result = process_line(line, line_num)
                if result:
                    outfile.write(result + '\n')
                    processed_count += 1
                else:
                    skipped_count += 1
                    # Error message already printed to stderr in process_line
                    pass
        
        # Print summary to stdout (normal output)
        print(f"Conversion completed:")
        print(f"  Successfully processed: {processed_count} lines")
        print(f"  Skipped: {skipped_count} lines")
        print(f"  Output file: {output_file}")
        # No separate error log file; all warnings go to stderr
        
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}", file=sys.stderr)
        sys.exit(1)
    except PermissionError as e:
        print(f"Error: Permission denied - {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: Unexpected error occurred - {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
