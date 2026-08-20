#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
Stat: basic statistics and structure plot for annotation results.
Supports -i (basic-info file), --rnaplot, and optional --list.
"""

import argparse
import os
import sys
import subprocess
from pathlib import Path
from datetime import datetime

def add_arguments(parser: argparse.ArgumentParser):
    """Define Stat-specific arguments."""
    parser.add_argument("-i", "--input", required=True,
                        help="Input basic-info file from annotation")
    parser.add_argument("-o", "--output",
                        help="Output directory (default: mirdeep-stat-<timestamp>)")
    parser.add_argument("--rnaplot", action="store_true",
                        help="Generate RNA secondary structure plots")
    parser.add_argument("--list",
                        help="Comma-separated list of miRNA IDs for RNA plotting (e.g., 'Sta-MIR157a,Sta-MIR156a')")
    parser.add_argument("-h", "--help", action="store_true",
                        help="Show this help message and exit.")


def run(args):
    """Execute Stat analysis."""
    # Resolve project root
    project_root = getattr(args, 'project_root', None) or Path(__file__).resolve().parents[2]
    scripts_dir = project_root / "scripts"

    # Validate input
    input_file = Path(args.input)
    if not input_file.is_file():
        sys.exit(f"Input basic-info file not found: {input_file}")

    # Output directory
    if args.output:
        output_dir = Path(args.output)
    else:
        output_dir = Path(f"mirdeep-stat-{datetime.now().strftime('%m%d%y-%H%M')}")
    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- Step 1: basic statistics ----
    basic_stat_script = scripts_dir / "basic_stat.R"
    cmd = f"Rscript {basic_stat_script} -i {input_file} -o {output_dir}"
    subprocess.run(cmd, shell=True, check=True)
    print("Basic statistics generated.")

    # ---- Step 2: RNA structure plot (if requested) ----
    if args.rnaplot:
        rna_plot_script = scripts_dir / "RNA_plot.py"
        if args.list:
            cmd = (f"python {rna_plot_script} -i {input_file} -o {output_dir} "
                   f"--list \"{args.list}\"")
        else:
            cmd = f"python {rna_plot_script} -i {input_file} -o {output_dir}"
        subprocess.run(cmd, shell=True, check=True)
        print("RNA structure plots generated.")

    print(f"\nStat analysis completed. Results in {output_dir}")
