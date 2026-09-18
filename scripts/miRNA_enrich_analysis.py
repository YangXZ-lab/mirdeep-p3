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

# Exit-code convention shared with enrich_analysis.R and the callers
# (functional_analysis.py / onestep.py):
#   0 success | 1 failure | 2 PARTIAL_SUCCESS (some sub-tasks failed)
EXIT_OK, EXIT_FAILURE, EXIT_PARTIAL = 0, 1, 2

# Exit code emitted by enrich_analysis.R when it ran fine but found no
# significant KEGG/GO term -- a legitimate empty result, not a failure.
R_NO_TERMS = 3


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

    # Process each miRNA, recording the outcome of every single task instead of
    # only printing the failure: a partly failed run must not exit 0.
    results = []          # (mirna, short, status, detail)

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
            rc = subprocess.run(cmd).returncode
        except OSError as e:
            # e.g. Rscript not on PATH -- record it, do not abort the batch
            rc = EXIT_FAILURE
            print(f"Error: could not run enrich_analysis.R for {mirna}: {e}",
                  file=sys.stderr)

        if rc == EXIT_OK:
            status, detail = "ok", f"results in {mirna_out}"
            print(f"Enrichment analysis finished for {mirna}; {detail}")
        elif rc == R_NO_TERMS:
            status = "no_terms"
            detail = "no significant terms for this gene list"
            print(f"Note: enrich_analysis.R found no significant terms for {mirna} "
                  f"(exit {rc}) -- empty result, not a failure")
        else:
            status, detail = "failed", f"enrich_analysis.R exit code {rc}"
            print(f"Error: enrich_analysis.R failed for {mirna} (code {rc})")

        results.append((mirna, short, status, detail))

    # Cleanup temporary files
    print("\nCleaning up temporary files...")
    try:
        shutil.rmtree(temp_dir)
        print(f"Temporary directory removed: {temp_dir}")
    except OSError as e:
        print(f"Warning: could not delete temp directory {temp_dir}: {e}",
              file=sys.stderr)

    # ---- Aggregate every sub-task exit code into one explicit verdict ----
    ok       = [r for r in results if r[2] == "ok"]
    no_terms = [r for r in results if r[2] == "no_terms"]
    failed   = [r for r in results if r[2] == "failed"]

    print("\n" + "=" * 72)
    print(f"miRNA enrichment summary: {len(results)} task(s) -- "
          f"{len(ok)} ok, {len(no_terms)} no significant terms, {len(failed)} failed")
    for mirna, short, status, detail in results:
        print(f"   {status.upper():<9} {mirna:<16} ({short})  {detail}")
    if no_terms:
        print("Note: 'no significant terms' is an empty result, not a failure -- it "
              "usually means too few of these genes carry GO/KEGG annotation; check "
              "the OrgDb mappings for these gene IDs before treating it as a bug.")
    print("=" * 72)

    if failed and len(failed) == len(results):
        print("STATUS: FAILURE -- every miRNA enrichment task failed; no results "
              "were produced.")
        return EXIT_FAILURE
    if failed:
        print(f"STATUS: PARTIAL_SUCCESS -- {len(failed)} of {len(results)} miRNA "
              f"task(s) failed; the other {len(ok) + len(no_terms)} finished. "
              f"Enrichment output is INCOMPLETE for: "
              f"{', '.join(r[0] for r in failed)}")
        return EXIT_PARTIAL
    print("STATUS: SUCCESS -- all miRNA enrichment tasks completed.")
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
