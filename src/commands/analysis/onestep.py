#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
Onestep: run a full analysis pipeline from basic-info to functional enrichment.
"""

import argparse
import os
import sys
import subprocess
import shutil
from pathlib import Path


def add_arguments(parser: argparse.ArgumentParser):
    parser.add_argument("-h", "--help", action="help",
                        help="Show this help message and exit.")
    parser.add_argument("-b", "--basic", required=True, help="Input basic-info file from annotation")
    parser.add_argument("-c", "--count", required=True, help="Input count file (mature.count)")
    parser.add_argument("-r", "--rpm", required=True, help="Input RPM file (mature.exp)")
    parser.add_argument("--fai", required=True, help="Genome FASTA index file (.fai)")
    parser.add_argument("-g", "--genome", required=True, help="Genome FASTA file")
    parser.add_argument("-t", "--transcript", required=True, help="Transcript FASTA file (for target prediction)")
    parser.add_argument("-p", "--protein", required=True, help="Protein FASTA file (for eggNOG)")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("--case1", required=True, help="Columns for group1 (comma-separated, e.g., 3,4,5)")
    parser.add_argument("--case2", required=True, help="Columns for group2 (comma-separated, e.g., 6,7,8)")
    parser.add_argument("--case1name", default="case1", help="Name of group1 (default: case1)")
    parser.add_argument("--case2name", default="case2", help="Name of group2 (default: case2)")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads (default: 1)")
    parser.add_argument("--rnaplot", action="store_true", help="Generate RNA structure plots in Stat")
    parser.add_argument("--tfbsplot", action="store_true", help="Generate TFBS report picture")
    parser.add_argument("--DEOnly", action="store_true", help="Focus only on DE miRNAs for downstream")
    parser.add_argument("--chord", action="store_true", help="Generate chord diagram in functional analysis")


def run(args):
    project_root = getattr(args, 'project_root', None) or Path(__file__).resolve().parents[2]
    main_script = project_root / "mirdeep-p3"
    python = sys.executable

    # ---- Validate case columns ----
    case1_cols = [int(x.strip()) for x in args.case1.split(",")]
    case2_cols = [int(x.strip()) for x in args.case2.split(",")]
    if not case1_cols or not case2_cols:
        sys.exit("Error: --case1 and --case2 must contain at least one column index each.")
    if any(c in (1, 2) for c in case1_cols + case2_cols):
        sys.exit("Error: column indices 1 or 2 are not allowed.")
    if set(case1_cols) & set(case2_cols):
        sys.exit("Error: --case1 and --case2 must not overlap.")
    if len(case1_cols) != len(case2_cols):
        sys.exit("Error: --case1 and --case2 must have the same number of columns.")
    
    output_dir = Path(args.output)
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- 1. Stat ----
    cmd = [python, str(main_script), "analysis", "Stat",
           "-i", args.basic, "-o", str(output_dir)]
    if args.rnaplot:
        cmd.append("--rnaplot")
    print("[onestep] Running Stat...")
    subprocess.run(cmd, check=True)

    # ---- 2. TFBS ----
    cmd = [python, str(main_script), "analysis", "TFBS",
           "-i", args.basic,
           "--fai", args.fai,
           "-g", args.genome,
           "-o", str(output_dir)]
    if args.tfbsplot:
        cmd.append("-p")
    print("[onestep] Running TFBS...")
    subprocess.run(cmd, check=True)

    # ---- 3. Target_finder ----
    cmd = [python, str(main_script), "analysis", "Target_finder",
           "-b", args.basic,
           "-c", args.transcript,
           "-o", str(output_dir),
           "-t", str(args.threads)]
    print("[onestep] Running Target_finder...")
    subprocess.run(cmd, check=True)

    # ---- 4. Differential_expression ----
    deg_cmd = [python, str(main_script), "analysis", "Differential_expression",
               "-c", args.count,
               "-r", args.rpm,
               "--case1", args.case1,
               "--case2", args.case2,
               "--case1name", args.case1name,
               "--case2name", args.case2name,
               "-o", str(output_dir)]
    if args.DEOnly:
        deg_cmd.append("--DEOnly")
    print("[onestep] Running Differential_expression...")
    subprocess.run(deg_cmd, check=True)

    # ---- 5. Functional_analysis ----
    if args.DEOnly:
        # Filter target_finder.tsv with DE miRNAs
        deg_res_file = output_dir / f"{args.case1name}_{args.case2name}_res.txt"
        if not deg_res_file.exists():
            sys.exit(f"DEG result file not found: {deg_res_file}")
        target_de = output_dir / "target_finder_DEG.tsv"
        awk_cmd = (
            f"cat {deg_res_file} | sed '1d' | awk '$NF!=\"NOT\"{{ print $1 }}' | "
            f"awk 'NR==FNR {{a[$1];next}} {{if($1 in a){{ print $0 }}}}' - {output_dir}/target_finder.tsv > {target_de}"
        )
        subprocess.run(awk_cmd, shell=True, check=True)
        target_input = str(target_de)
    else:
        target_input = str(output_dir / "target_finder.tsv")

    func_cmd = [python, str(main_script), "analysis", "Functional_analysis",
                "-p", args.protein,
                "--target", target_input,
                "-t", str(args.threads),
                "-o", str(output_dir)]
    if args.chord:
        func_cmd.append("--chord")
    print("[onestep] Running Functional_analysis...")
    subprocess.run(func_cmd, check=True)

    # ---- Cleanup temp ----
    temp_dir = output_dir / "temp"
    if temp_dir.exists():
        shutil.rmtree(temp_dir)

    print("Onestep pipeline completed successfully.")