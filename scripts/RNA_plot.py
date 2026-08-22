#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
Convert TSV of RNA sequences and structure annotations to SVG plots via RNAplot.
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate RNA secondary structure SVGs from a TSV file."
    )
    parser.add_argument("-i", "--input", required=True, help="Input TSV file")
    parser.add_argument(
        "-o", "--output", default="output", help="Output directory (default: output)"
    )
    parser.add_argument(
        "--list",
        default=None,
        help="Comma-separated miRNA names to process. If omitted, process all.",
    )
    return parser.parse_args()


def process_row(row, offset, output_dir):
    """Process a single TSV row: write temp struct, call RNAplot, return True on success."""
    name = row[0].strip()
    seq = row[8]
    stru = row[9]

    # Convert absolute genomic coordinates to 1‑based positions relative to offset ($7)
    rel = lambda x: int(x) - offset + 1
    mat_s  = rel(row[15])   # column $16
    mat_e  = rel(row[16])   # column $17
    star_s = rel(row[19])   # column $20
    star_e = rel(row[20])   # column $21

    # Write temporary struct file
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".struct", delete=False, encoding="utf-8"
    ) as tmpf:
        tmpf.write(f">{name}\n{seq}\n{stru}\n")
        tmp_path = tmpf.name

    try:
        # Correct --pre format: start end color [type]
        # NOTE: 'Fomark' is not a valid type; use default rectangle (or 'rect').
        annotation = f"{mat_s} {mat_e} red {star_s} {star_e} blue"

        with open(tmp_path) as infile:
            subprocess.run(
                [
                    "RNAplot",
                    "--output-format=svg",
                    "--pre",
                    annotation,
                ],
                stdin=infile,
                cwd=output_dir,
                check=True,
                capture_output=True,
                text=True,
            )
        print(f"✅ Done: {name} → {output_dir}/{name}_ss.svg")
        return True

    except subprocess.CalledProcessError as e:
        print(f"❌ RNAplot failed for {name}: {e.stderr.strip()}", file=sys.stderr)
        return False

    finally:
        os.unlink(tmp_path)


def main():
    args = parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build optional name filter
    name_filter = None
    if args.list is not None:
        name_filter = {n.strip() for n in args.list.split(",") if n.strip()}
        if not name_filter:
            print("Warning: --list provided but empty. No output will be produced.")

    try:
        with open(args.input, newline="", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                # Skip empty lines
                if not row or all(cell.strip() == "" for cell in row):
                    continue
                # Need at least 21 columns (indices 0..20)
                if len(row) < 21:
                    print(
                        f"Warning: row has {len(row)} columns, need ≥21. Skipped: {row}",
                        file=sys.stderr,
                    )
                    continue

                name = row[0].strip()
                if name_filter and name not in name_filter:
                    continue

                # offset = column $7 (0‑based index 6)
                try:
                    offset = int(row[6])
                except ValueError:
                    print(
                        f"Warning: invalid offset value for {name}: column 7 = {row[6]}. Skipped.",
                        file=sys.stderr,
                    )
                    continue

                process_row(row, offset, str(output_dir))

    except FileNotFoundError:
        sys.exit(f"Error: input file not found: {args.input}")
    except Exception as e:
        sys.exit(f"Unexpected error: {e}")


if __name__ == "__main__":
    main()
