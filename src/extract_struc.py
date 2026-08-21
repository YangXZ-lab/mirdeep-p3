#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
extract_struc.py

Predict RNA secondary structures for miRDP3 filtered annotations.

Usage:
    python extract_struc.py -i filtered-all-nr-anno \\
        -genome /path/to/genome.fa \\
        --threads 40 \\
        -struc filtered-all-nr-anno-pre.struc \\
        -struc_20nt stemloop_20nt.struc \\
        -fasta_20nt stemloop_20nt.fasta
"""

import sys
import argparse
import subprocess
import tempfile
import os

def run_rnafold_precursor(input_file, output_file, threads):
    """
    Step 1: Extract col4 and col8, run RNAfold.
    """
    # Use awk to format: ">" + col4 + "\n" + col8
    awk_cmd = ["awk", '{print ">"$4"\\n"$8}', input_file]
    rnafold_cmd = ["RNAfold", "--noPS", f"-j{threads}"]

    with open(output_file, 'w') as out_fh:
        p1 = subprocess.Popen(awk_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2 = subprocess.Popen(rnafold_cmd, stdin=p1.stdout, stdout=out_fh, stderr=subprocess.PIPE)
        p1.stdout.close()
        _, stderr = p2.communicate()
        if p2.returncode != 0:
            sys.stderr.write(f"Error in RNAfold (precursor): {stderr.decode()}\n")
            sys.exit(1)

def generate_bed_and_fetch(input_file, genome_file, output_fasta, threads):
    """
    Step 2: Generate BED intervals extended by 20 nt and run bedtools getfasta.
    """
    # Step 2a: Extract fields and replace ".." with tab
    awk_extract = ["awk", '{print $1"\t"$6"\t"$4"\t"$3"\t"$2}', input_file]
    sed_cmd = ["sed", 's/\\.\\./\t/g']
    # Step 2b: Compute extended coordinates (start-20, end+20) and output 6 columns
    awk_extend = ["awk", '{print $1"\t"($2-20)"\t"($3+20)"\t"$4"\t"$5"\t"$6}']

    bedtools_cmd = [
        "bedtools", "getfasta",
        "-s", "-nameOnly",
        "-bed", "-",
        "-fi", genome_file,
        "-fo", output_fasta
    ]

    p1 = subprocess.Popen(awk_extract, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen(sed_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p1.stdout.close()
    p3 = subprocess.Popen(awk_extend, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2.stdout.close()
    p4 = subprocess.Popen(bedtools_cmd, stdin=p3.stdout, stderr=subprocess.PIPE)
    p3.stdout.close()

    _, stderr = p4.communicate()
    if p4.returncode != 0:
        sys.stderr.write(f"Error in bedtools getfasta: {stderr.decode()}\n")
        sys.exit(1)

def run_rnafold_20nt(input_fasta, output_file, threads):
    """
    Step 3: Run RNAfold on the extracted 20nt-extended sequences.
    """
    rnafold_cmd = ["RNAfold", "--noPS", f"-j{threads}", input_fasta]
    with open(output_file, 'w') as out_fh:
        proc = subprocess.run(rnafold_cmd, stdout=out_fh, stderr=subprocess.PIPE)
        if proc.returncode != 0:
            sys.stderr.write(f"Error in RNAfold (20nt): {proc.stderr.decode()}\n")
            sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Predict RNA secondary structures for miRDP2 filtered annotations."
    )
    parser.add_argument('-i', '--input', required=True,
                        help='Input filtered annotation file (filtered-all-nr-anno)')
    parser.add_argument('-genome', '--genome', required=True,
                        help='Genome FASTA file for sequence extraction')
    parser.add_argument('--threads', type=int, required=True,
                        help='Number of threads for RNAfold (used as -j)')
    parser.add_argument('-struc', '--struc_out', required=True,
                        help='Output RNAfold result for precursor sequences')
    parser.add_argument('-struc_20nt', '--struc_20nt_out', required=True,
                        help='Output RNAfold result for 20nt-extended sequences')
    parser.add_argument('-fasta_20nt', '--fasta_20nt_out', required=True,
                        help='Output FASTA file for 20nt-extended sequences')
    args = parser.parse_args()

    # Check that required tools are available
    try:
        subprocess.run(["RNAfold", "--version"], capture_output=True, check=True)
    except (subprocess.SubprocessError, FileNotFoundError):
        sys.stderr.write("Error: RNAfold not found in PATH.\n")
        sys.exit(1)

    try:
        subprocess.run(["bedtools", "--version"], capture_output=True, check=True)
    except (subprocess.SubprocessError, FileNotFoundError):
        sys.stderr.write("Error: bedtools not found in PATH.\n")
        sys.exit(1)

    # Step 1
    print("Step 1: Running RNAfold on precursor sequences...", file=sys.stderr)
    run_rnafold_precursor(args.input, args.struc_out, args.threads)

    # Step 2
    print("Step 2: Extracting 20nt-extended genomic intervals...", file=sys.stderr)
    generate_bed_and_fetch(args.input, args.genome, args.fasta_20nt_out, args.threads)

    # Step 3
    print("Step 3: Running RNAfold on 20nt-extended sequences...", file=sys.stderr)
    run_rnafold_20nt(args.fasta_20nt_out, args.struc_20nt_out, args.threads)

    print("All steps completed successfully.", file=sys.stderr)

if __name__ == "__main__":
    main()
