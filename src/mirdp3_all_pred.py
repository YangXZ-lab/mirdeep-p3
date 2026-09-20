#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
Convert P_prediction files to RPM-normalized format with modified identifiers.

"""

import sys
import os
import re
import argparse
from typing import Dict, List, TextIO, Optional


def extract_sample_name(pred_file: str) -> Optional[str]:
    """
    Extract sample name from a prediction file path.
    
    The file name should match either of the patterns:
        * <sample>_filter_P_prediction
        * <sample>_P_prediction
    
    Args:
        pred_file: Path to the prediction file.
    
    Returns:
        Sample name if matched, otherwise None.
    """
    base = os.path.basename(pred_file)
    # Try pattern with '_filter_P_prediction' first (original Perl regex)
    match = re.match(r'^(\S+)_filter_P_prediction$', base)
    if match:
        return match.group(1)
    # Fallback to pattern without 'filter'
    match = re.match(r'^(\S+)_P_prediction$', base)
    if match:
        return match.group(1)
    return None


def read_counts_file(counts_file: str) -> Dict[str, int]:
    """
    Read the total read counts per sample.
    
    Expected format: tab-separated, two columns: sample_name<TAB>total_reads
    
    Args:
        counts_file: Path to the counts file.
    
    Returns:
        Dictionary mapping sample name -> total reads (int).
    
    Raises:
        SystemExit: If file cannot be read or format is invalid.
    """
    counts = {}
    try:
        with open(counts_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    print(f"Warning: Line {line_num} in {counts_file} "
                          f"has <2 columns, skipping: {line}", file=sys.stderr)
                    continue
                sample, total_str = parts[0], parts[1]
                try:
                    counts[sample] = int(total_str)
                except ValueError:
                    print(f"Warning: Invalid total reads '{total_str}' for "
                          f"sample '{sample}' at line {line_num}, skipping.",
                          file=sys.stderr)
    except Exception as e:
        print(f"Error: Cannot read counts file {counts_file}: {e}",
              file=sys.stderr)
        sys.exit(1)
    
    if not counts:
        print(f"Error: No valid sample counts found in {counts_file}",
              file=sys.stderr)
        sys.exit(1)
    
    return counts


def process_prediction_file(pred_file: str,
                            sample_name: str,
                            total_reads: int,
                            output_handle: TextIO,
                            map_handle: Optional[TextIO],
                            series_start: int = 1) -> int:
    """
    Process a single prediction file, apply transformations, and write to output.
    
    Transformations per line:
        1. Replace the 4th column (index 3) with:
           "{sample_name}-{original_4th_col}-{3rd_col}-{line_series}"
        2. Extract the integer after '_x' in the 3rd column as raw read count.
        3. Compute RPM = (raw_reads / total_reads) * 1,000,000.
        4. Replace every occurrence of '_x{raw_reads}' with '_x{RPM}'.
    
    Args:
        pred_file: Path to the input prediction file.
        sample_name: Name of the sample (used in identifier prefix).
        total_reads: Total number of reads for this sample.
        output_handle: Open file handle to write transformed lines.
        map_handle: Open file handle to write mapping (original_x -> rpm_x), or None.
        series_start: Starting series number (usually 1 per file).
    
    Returns:
        The last series number used (series_start + lines_processed - 1).
    """
    series = series_start
    try:
        with open(pred_file, 'r') as f:
            for line in f:
                line = line.rstrip('\n')
                if not line:
                    continue  # skip empty lines
                
                # Split into columns (tab-separated)
                cols = line.split('\t')
                if len(cols) < 4:
                    print(f"Warning: Skipping malformed line in {pred_file}: "
                          f"expected at least 4 columns, got {len(cols)}",
                          file=sys.stderr)
                    continue
                
                # Original values needed for transformation
                col_2_original = cols[2]    # e.g., read01402880_x1045
                col_3 = cols[3]             # e.g., Chr05_164
                
                # --- Step 1: Replace 4th column with new identifier ---
                new_identifier = f"{sample_name}-{col_3}-{col_2_original}-{series}"
                escaped_col3 = re.escape(col_3)
                line = re.sub(escaped_col3, new_identifier, line, count=1)
                
                # --- Step 2: Extract raw read count from col_2 ---
                read_match = re.search(r'_x(\d+)', col_2_original)
                if not read_match:
                    print(f"Warning: Cannot extract read count from '{col_2_original}' "
                          f"in {pred_file}, skipping line", file=sys.stderr)
                    series += 1
                    continue
                raw_reads = int(read_match.group(1))
                
                # --- Step 3: Calculate RPM ---
                rpm = (raw_reads / total_reads) * 1_000_000
                rpm_str = str(rpm)
                
                # --- Step 4: Replace all '_x{raw_reads}' with '_x{rpm_str}' ---
                line = line.replace(f"_x{raw_reads}", f"_x{rpm_str}")
                
                # Write transformed line
                output_handle.write(line + '\n')
                
                # --- Optional: Write mapping from original read ID to RPM version ---
                if map_handle is not None:
                    # Extract the new read ID from the transformed line (second column)
                    new_cols = line.split('\t')
                    if len(new_cols) >= 3:
                        col_2_new = new_cols[2]
                        map_handle.write(f"{col_2_original}\t{col_2_new}\n")
                    else:
                        print(f"Warning: Cannot extract new read ID from transformed line",
                              file=sys.stderr)
                
                series += 1
                
    except Exception as e:
        print(f"Error processing file {pred_file}: {e}", file=sys.stderr)
        # Continue with next file, but stop processing this one
    return series - 1


def main():
    """Command line interface."""
    parser = argparse.ArgumentParser(
        description='Convert P_prediction files to RPM-normalized format.',
        epilog='Example: %(prog)s -i sample1_P_prediction,sample2_P_prediction '
               '-c counts.txt -o all_predictions.rpm.txt -m mapping.txt',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Comma-separated list of P_prediction files to process.'
    )
    
    parser.add_argument(
        '-c', '--counts',
        required=True,
        help='File containing total reads per sample (two columns: sample<TAB>total_reads).'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output file path. Transformed lines from all input files will be '
             'written here (overwrites by default, use -a to append).'
    )
    
    parser.add_argument(
        '-m', '--map',
        help='Optional mapping file: each line contains original read ID (with _xCOUNT) '
             'and converted read ID (with _xRPM), tab-separated.'
    )
    
    parser.add_argument(
        '-a', '--append',
        action='store_true',
        help='Append to output file instead of overwriting.'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print detailed progress information.'
    )
    
    args = parser.parse_args()
    
    # Parse input file list
    input_files = [f.strip() for f in args.input.split(',') if f.strip()]
    if not input_files:
        print("Error: No input files specified.", file=sys.stderr)
        sys.exit(1)
    
    # Read total counts per sample
    counts = read_counts_file(args.counts)
    if args.verbose:
        print(f"Read counts for {len(counts)} sample(s):", file=sys.stderr)
        for samp, cnt in counts.items():
            print(f"  {samp}: {cnt:,}", file=sys.stderr)
    
    # Open output file (append or write)
    mode = 'a' if args.append else 'w'
    try:
        out_fh = open(args.output, mode)
    except Exception as e:
        print(f"Error: Cannot open output file {args.output}: {e}",
              file=sys.stderr)
        sys.exit(1)
    
    # Open mapping file if requested (always overwrite, as it's a derived mapping)
    map_fh = None
    if args.map:
        try:
            map_fh = open(args.map, 'w')
            # Write header (optional)
            map_fh.write("#original_read_id\trpm_converted_read_id\n")
        except Exception as e:
            print(f"Error: Cannot open mapping file {args.map}: {e}",
                  file=sys.stderr)
            sys.exit(1)
    
    processed_count = 0
    skipped_count = 0
    
    with out_fh:
        # If map_fh is not None, it will be closed automatically after the block
        with map_fh if map_fh else open(os.devnull, 'w') as mfh:
            for pred_file in input_files:
                if not os.path.isfile(pred_file):
                    print(f"Warning: Input file not found, skipping: {pred_file}",
                          file=sys.stderr)
                    skipped_count += 1
                    continue
                
                # Extract sample name from filename
                sample_name = extract_sample_name(pred_file)
                if sample_name is None:
                    print(f"Warning: Cannot extract sample name from {pred_file}, "
                          f"skipping (expected pattern: *_P_prediction or *_filter_P_prediction).",
                          file=sys.stderr)
                    skipped_count += 1
                    continue
                
                # Look up total reads for this sample
                total_reads = counts.get(sample_name)
                if total_reads is None:
                    print(f"Warning: Sample '{sample_name}' not found in counts file, "
                          f"skipping {pred_file}.", file=sys.stderr)
                    skipped_count += 1
                    continue
                
                if args.verbose:
                    print(f"Processing {pred_file} (sample: {sample_name}, "
                          f"total reads: {total_reads:,})", file=sys.stderr)
                
                # Process the file
                last_series = process_prediction_file(
                    pred_file, sample_name, total_reads, out_fh, map_fh, series_start=1
                )
                processed_count += 1
                if args.verbose:
                    print(f"  -> Processed {last_series} lines.", file=sys.stderr)
    
    # Summary
    print(f"\nProcessing completed.", file=sys.stderr)
    print(f"  Processed files: {processed_count}", file=sys.stderr)
    if skipped_count:
        print(f"  Skipped files: {skipped_count}", file=sys.stderr)
    print(f"  Output written to: {args.output}", file=sys.stderr)
    if args.map:
        print(f"  Mapping written to: {args.map}", file=sys.stderr)


if __name__ == "__main__":
    main()
