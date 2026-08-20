#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
Target_finder: miRNA target prediction (ssearch36 based).
Supports input via -i (FASTA) or -b (basic-info file from annotation).
"""

import argparse
import os
import sys
import subprocess
import shutil
from pathlib import Path
from datetime import datetime

from utils.dependencies import check_external_tools


def add_arguments(parser: argparse.ArgumentParser):
    """Define Target_finder-specific arguments."""
    parser.add_argument("-i", "--input", help="Input mature miRNA FASTA file")
    parser.add_argument("-b", "--basic", help="Input basic-info file from annotation (instead of -i)")
    parser.add_argument("-c", "--cds", required=True, help="CDS sequences FASTA file")
    parser.add_argument("-o", "--output", help="Output directory / prefix (def. mirdeep-target_finder-<timestamp>)")
    parser.add_argument("-e", "--evalue", type=float, default=2.5,
                        help="E-value cutoff for target filtering (default: 2.5)")
    parser.add_argument("--GUs", type=int, default=1,
                        help="Maximum number of G:U wobbles (default: 1)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads for ssearch36 (default: 1)")


def run(args):
    """Execute Target_finder analysis."""
    # ---- Resolve project root ----
    project_root = getattr(args, 'project_root', None) or Path(__file__).resolve().parents[2]
    src_dir = project_root / "src"
    mi_rna_target_dir = src_dir / "MiRNATarget"
    scripts_dir = project_root / "scripts"

    # ---- Determine input ----
    if args.input and args.basic:
        sys.exit("Error: please specify either -i (FASTA) or -b (basic-info), not both.")
    if not args.input and not args.basic:
        sys.exit("Error: you must provide either -i or -b.")

    # ---- Output directory ----
    if args.output:
        output_dir = Path(args.output)
    else:
        output_dir = Path(f"mirdeep-target_finder-{datetime.now().strftime('%m%d%y-%H%M')}")
    output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = output_dir / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)

    # ---- Prepare mature FASTA ----
    if args.basic:
        basic_info_file = Path(args.basic).resolve()
        if not basic_info_file.is_file():
            sys.exit(f"Basic-info file not found: {basic_info_file}")
        mature_fa = temp_dir / "mature.fa"
        # Extract mature sequences: column 1 (ID) and column 18 (sequence)
        cmd = f"cat {basic_info_file} | awk '{{ print \">\"$1\"\\n\"$18 }}' > {mature_fa}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"Extracted mature sequences from {basic_info_file} into {mature_fa}")
    else:
        mature_fa = Path(args.input).resolve()
        if not mature_fa.is_file():
            sys.exit(f"Input FASTA file not found: {mature_fa}")

    # ---- Check required external tools ----
    tools_ok, missing = check_external_tools({'ssearch36': None})
    if not tools_ok:
        sys.exit("Error: ssearch36 not found in PATH.")

    # ---- Capture common parameters ----
    threads = args.threads

    # ---- Step 1: Reverse complement ----
    rp_fa = temp_dir / "mature.rp.fa"
    subprocess.run(f"seqkit seq {mature_fa} -r -p > {rp_fa}", shell=True, check=True)

    # ---- Step 2: Forward and reverse alignments ----
    cds_fa = Path(args.cds).resolve()
    if not cds_fa.exists():
        sys.exit(f"CDS file not found: {cds_fa}")

    # ssearch36 base options (threads now comes from -t)
    ssearch_opts = (f"-f -8 -g -3 -E 10000 -T {threads} -b 200 "
                    "-r +4/-3 -n -U -W 10 -N 20000")

    # Forward alignment
    forw_aln = temp_dir / "forward_alignment"
    subprocess.run(f"ssearch36 {ssearch_opts} {mature_fa} {cds_fa} > {forw_aln}",
                   shell=True, check=True)

    # Reverse alignment
    rev_aln = temp_dir / "reverse_alignment"
    subprocess.run(f"ssearch36 {ssearch_opts} {rp_fa} {cds_fa} > {rev_aln}",
                   shell=True, check=True)

    # ---- Step 3: Parse alignments to TSV ----
    parse_script = scripts_dir / "parse_ssearch.py"
    forw_tsv = temp_dir / "forward_alignment.tsv"
    rev_tsv = temp_dir / "reverse_alignment.tsv"
    subprocess.run(f"python {parse_script} -i {forw_aln} > {forw_tsv}", shell=True, check=True)
    subprocess.run(f"python {parse_script} -i {rev_aln} > {rev_tsv}", shell=True, check=True)
    total_tsv = temp_dir / "total.tsv"
    subprocess.run(f"cat {forw_tsv} {rev_tsv} > {total_tsv}", shell=True, check=True)

    # ---- Step 4: Statistical filtering ----
    parse_targets_script = mi_rna_target_dir / "parse_mirna_targets.py"
    raw_target = temp_dir / "target_finder_raw.tsv"
    cmd = (f"python {parse_targets_script} -i {total_tsv} "
           f"--E_cutoff {args.evalue} --GUs_cutoff {args.GUs} > {raw_target}")
    subprocess.run(cmd, shell=True, check=True)

    # Get mature length
    mature_len = temp_dir / "mature.length"
    subprocess.run(f"seqkit fx2tab --length --name --header-line {mature_fa} | grep -v '#' > {mature_len}",
                   shell=True, check=True)

    # Final filtered output
    final_output = output_dir / "target_finder.tsv"
    subprocess.run(
        f"cat {raw_target} | grep -v '#' | "
        f"awk 'NR==FNR {{a[$1]=$2;next}} {{if($1 in a){{ print $0\"\\t\"a[$1] }}}}' {mature_len} - | "
        f"awk '{{if($(NF-6)>=($NF-1)){{ print $0 }}}}' | "
        f"awk -F \"\\t\" '{{ print $1\"\\t\"$2\"\\t\"$NF\"\\t\"$8\"\\t\"($8+$NF-1)\"\\t\"$(NF-2) }}' > {final_output}",
        shell=True, check=True)

    # ---- Step 5: Clean temp ----
    # shutil.rmtree(temp_dir)

    print(f"\nTarget_finder completed. Results written to {final_output}")
