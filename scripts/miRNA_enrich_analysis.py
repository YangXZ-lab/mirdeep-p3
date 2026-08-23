#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
Pipeline to orchestrate miRNA target gene enrichment analysis.

This script:
  1. Reads an input file containing miRNA-target gene relationships.
  2. Extracts unique miRNAs and their associated target genes.
  3. Creates an output subdirectory for each miRNA (using the suffix after '-').
  4. Writes a temporary gene list for each miRNA and invokes the companion
     R script (enrich_analysis.R) for enrichment analysis.
  5. Cleans up temporary files automatically.

Usage example:
  python mirna_enrich_pipeline.py \
      -i input.tsv \
      --orgdb /path/to/orgdb/org.Morg.eg.db \
      -f /path/to/orgdb_dir \
      -o /path/to/output
"""

import argparse
import os
import sys
import subprocess
import shutil


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Orchestrate miRNA target enrichment analysis"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input file; first two columns are miRNA and target gene ID"
    )
    parser.add_argument(
        "--orgdb",
        required=True,
        help="Path to the OrgDB database, passed to enrich_analysis.R via -i"
    )
    parser.add_argument(
        "-f", "--file",
        required=True,
        help="Directory containing auxiliary files, passed to enrich_analysis.R via -f"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Root output directory for results and temporary files"
    )
    return parser.parse_args()


def read_mirna_genes(input_path):
    """
    Build a dict mapping miRNA to a set of target genes from the input file.

    Args:
        input_path (str): Path to the input file.

    Returns:
        dict: {mirna_name: set(genes)}.
    """
    mirna_genes = {}
    with open(input_path, 'r') as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                print(f"Warning: skipping line {lineno}, not enough columns: {line}",
                      file=sys.stderr)
                continue
            mirna, gene = parts[0], parts[1]
            mirna_genes.setdefault(mirna, set()).add(gene)
    return mirna_genes


def get_short_name(full_name):
    """
    Extract the short miRNA name used for directory creation.

    The short name is the substring after the first '-'.
    If no '-' is present, the full name is returned.
    """
    return full_name.split('-', 1)[1] if '-' in full_name else full_name


def main():
    args = parse_args()

    # Validate input file
    if not os.path.isfile(args.input):
        sys.exit(f"Error: input file not found: {args.input}")

    # Locate the companion R script (same directory as this script)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script = os.path.join(script_dir, "enrich_analysis.R")
    if not os.path.isfile(r_script):
        sys.exit(f"Error: enrich_analysis.R not found at {r_script}")

    # Parse miRNA-gene relationships
    print(f"Reading input: {args.input}")
    mirna_genes = read_mirna_genes(args.input)
    if not mirna_genes:
        sys.exit("No miRNA-target gene relationships found. Exiting.")

    unique_mirnas = sorted(mirna_genes.keys())
    print(f"Unique miRNAs detected ({len(unique_mirnas)}): {', '.join(unique_mirnas)}")

    # Prepare output infrastructure
    os.makedirs(args.output, exist_ok=True)
    temp_dir = os.path.join(args.output, "temp")
    os.makedirs(temp_dir, exist_ok=True)

    # Process each miRNA
    for mirna, genes in mirna_genes.items():
        short = get_short_name(mirna)
        print(f"\n>>> Processing {mirna} (short: {short})")

        mirna_out = os.path.join(args.output, short)
        os.makedirs(mirna_out, exist_ok=True)

        # Write temporary gene list
        gene_file = os.path.join(temp_dir, f"{short}_genes.txt")
        with open(gene_file, 'w') as fh:
            for gene in sorted(genes):
                fh.write(f"{gene}\n")
        print(f"Temporary gene list written: {gene_file} ({len(genes)} genes)")

        # Build R command
        cmd = [
            "Rscript", r_script,
            "-i", args.orgdb,
            "-f", args.file,
            "-g", gene_file,
            "-o", mirna_out
        ]
        print(f"Running: {' '.join(cmd)}")

        try:
            subprocess.run(cmd, check=True)
            print(f"Enrichment analysis finished for {mirna}; results in {mirna_out}")
        except subprocess.CalledProcessError as e:
            print(f"Error: enrich_analysis.R failed for {mirna} (code {e.returncode})")
            # Continue with remaining miRNAs

    # Cleanup temporary files
    print("\nCleaning up temporary files...")
    try:
        shutil.rmtree(temp_dir)
        print(f"Temporary directory removed: {temp_dir}")
    except OSError as e:
        print(f"Warning: could not delete temp directory {temp_dir}: {e}",
              file=sys.stderr)

    print("All miRNA enrichment tasks completed.")


if __name__ == "__main__":
    main()
