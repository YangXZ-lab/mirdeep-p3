#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
Differential_expression: miRNA differential expression analysis.
Supports PCA, correlation, DEG, and optional dotplots.
Now includes an expression filter: miRNAs with expression < 5 in all samples
are removed before analysis.
"""

import argparse
import os
import sys
import subprocess
import shutil
from pathlib import Path
from datetime import datetime


def add_arguments(parser: argparse.ArgumentParser):
    """Define Differential_expression-specific arguments."""
    parser.add_argument("-h", "--help", action="help",
                        help="Show this help message and exit.")
    parser.add_argument("-c", "--count", required=True,
                        help="Input count matrix file (e.g., mature.count)")
    parser.add_argument("-r", "--rpm", required=True,
                        help="Input RPM/expression matrix file (e.g., mature.exp)")
    parser.add_argument("--case1", required=True,
                        help="Column indices for group1 (comma-separated, e.g., 3,4,5)")
    parser.add_argument("--case2", required=True,
                        help="Column indices for group2 (comma-separated, e.g., 6,7,8)")
    parser.add_argument("--case1name", default="case1",
                        help="Name of group1 (default: case1)")
    parser.add_argument("--case2name", default="case2",
                        help="Name of group2 (default: case2)")
    parser.add_argument("-o", "--output",
                        help="Output directory (default: mirdeep-de-<timestamp>)")
    parser.add_argument("--miRNA",
                        help="Comma-separated list of miRNA names to plot (optional)")
    parser.add_argument("-f", "--file",
                        help="File containing miRNA names, one per line (optional)")
    parser.add_argument("--DEOnly", action="store_true",
                        help="Plot only differentially expressed miRNAs")
    parser.add_argument("--min-expr", type=float, default=5.0,
                        help="Minimum expression in at least one sample to retain miRNA (default: 5.0)")


def filter_low_expression(exp_file, count_file, threshold=5.0):
    """
    Remove miRNAs that do not reach the given expression threshold in at least
    one sample. Both the expression and the count file are overwritten in place
    with only the retained miRNAs.
    """
    # Convert pathlib.Path to str if necessary
    exp_file = str(exp_file)
    count_file = str(count_file)

    # Read expression file and determine which miRNAs to keep
    kept_mirnas = []
    with open(exp_file, 'r') as f:
        header = f.readline().rstrip('\n')
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            mirna = parts[0]
            # Check if any value >= threshold
            values = []
            for v in parts[1:]:
                try:
                    values.append(float(v))
                except ValueError:
                    pass
            if any(v >= threshold for v in values):
                kept_mirnas.append(mirna)

    if not kept_mirnas:
        sys.exit("Error: no miRNA passes the expression filter.")

    kept_set = set(kept_mirnas)

    # Filter expression file
    tmp_exp = exp_file + ".tmp"
    with open(exp_file, 'r') as fin, open(tmp_exp, 'w') as fout:
        header = fin.readline()
        fout.write(header)
        for line in fin:
            if not line.strip():
                continue
            mir = line.split('\t', 1)[0]
            if mir in kept_set:
                fout.write(line)
    shutil.move(tmp_exp, exp_file)

    # Filter count file
    tmp_cnt = count_file + ".tmp"
    with open(count_file, 'r') as fin, open(tmp_cnt, 'w') as fout:
        header = fin.readline()
        fout.write(header)
        for line in fin:
            if not line.strip():
                continue
            mir = line.split('\t', 1)[0]
            if mir in kept_set:
                fout.write(line)
    shutil.move(tmp_cnt, count_file)

    return len(kept_mirnas)


def run(args):
    """Execute differential expression analysis."""
    project_root = getattr(args, 'project_root', None) or Path(__file__).resolve().parents[2]
    scripts_dir = project_root / "scripts"

    # ---- 1. Validate parameters ----
    case1_cols = [int(x.strip()) for x in args.case1.split(",")]
    case2_cols = [int(x.strip()) for x in args.case2.split(",")]

    if not case1_cols or not case2_cols:
        sys.exit("Error: --case1 and --case2 must contain at least one column index each.")

    # Check the first file's header to decide if column 2 is forbidden
    rpm_file = Path(args.rpm).resolve()
    if not rpm_file.is_file():
        sys.exit(f"RPM file not found: {rpm_file}")
    with open(rpm_file) as f:
        header_line = f.readline().strip()
    header_cols = header_line.split('\t')

    # Always forbid column 1 (miRNA names)
    forbidden_cols = {1}
    # If the second column is exactly "MIR_ID", also forbid column 2
    if len(header_cols) > 1 and header_cols[1] == "MIR_ID":
        forbidden_cols.add(2)

    all_cols = case1_cols + case2_cols
    bad = [c for c in all_cols if c in forbidden_cols]
    if bad:
        sys.exit(f"Error: column indices {', '.join(str(c) for c in bad)} are not allowed. "
                 f"Column 1 is always forbidden. Column 2 is forbidden only if the header is 'MIR_ID'.")

    if set(case1_cols) & set(case2_cols):
        sys.exit("Error: --case1 and --case2 must not overlap.")
    if len(case1_cols) != len(case2_cols):
        sys.exit("Error: --case1 and --case2 must have the same number of columns (replicates).")

    replicates = len(case1_cols)

    # miRNA list conflicts
    if args.miRNA and args.file:
        sys.exit("Error: --miRNA and -f/--file cannot be used together.")
    if (args.miRNA or args.file) and args.DEOnly:
        sys.exit("Error: --DEOnly cannot be combined with --miRNA or -f/--file.")

    # ---- 2. Output directory ----
    if args.output:
        output_dir = Path(args.output)
    else:
        output_dir = Path(f"mirdeep-de-{datetime.now().strftime('%m%d%y-%H%M')}")
    output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = output_dir / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)

    # ---- 3. Validate input files ----
    count_file = Path(args.count).resolve()
    if not count_file.is_file():
        sys.exit(f"Count file not found: {count_file}")

    # ---- 4. Extract columns ----
    extract_script = scripts_dir / "extract_col.py"

    # Process RPM
    tmp_exp = temp_dir / f"{rpm_file.stem}-temp{rpm_file.suffix}"
    cmd = (f"python {extract_script} -i {rpm_file} "
           f"--case1 {args.case1} --case2 {args.case2} -o {tmp_exp}")
    subprocess.run(cmd, shell=True, check=True)

    # Process count
    tmp_count = temp_dir / f"{count_file.stem}-temp{count_file.suffix}"
    cmd = (f"python {extract_script} -i {count_file} "
           f"--case1 {args.case1} --case2 {args.case2} -o {tmp_count}")
    subprocess.run(cmd, shell=True, check=True)

    # ---- 4b. Expression filter ----
    retained = filter_low_expression(tmp_exp, tmp_count, threshold=args.min_expr)
    print(f"Retained {retained} miRNAs after expression filter (>= {args.min_expr} in at least one sample).")

    # ---- 5. PCA ----
    pca_script = scripts_dir / "pca_miRNA.R"
    cmd = f"Rscript {pca_script} -i {tmp_exp} -r {replicates} -o {output_dir}"
    subprocess.run(cmd, shell=True, check=True)
    print("PCA plot generated.")

    # ---- 6. Correlation ----
    cor_script = scripts_dir / "correlation_miRNA.R"
    cmd = f"Rscript {cor_script} -i {tmp_exp} -o {output_dir}"
    subprocess.run(cmd, shell=True, check=True)
    print("Correlation plot generated.")

    # ---- 7. Differential expression ----
    deg_script = scripts_dir / "miRNA_DEG.R"
    cmd = (f"Rscript {deg_script} -i {tmp_count} -r {replicates} "
           f"-o {output_dir} --case1 \"{args.case1name}\" --case2 \"{args.case2name}\"")
    subprocess.run(cmd, shell=True, check=True)
    print("Differential expression analysis completed.")

    # ---- 8. miRNA dotplot ----
    dotplot_script = scripts_dir / "miRNA_dotplot.R"

    if args.DEOnly:
        deg_result_file = output_dir / f"{args.case1name}_{args.case2name}_res.txt"
        if not deg_result_file.exists():
            sys.exit(f"DEG result file not found: {deg_result_file}")
        de_mirna_list = temp_dir / "miRNA.list"
        subprocess.run(
            f"cat {deg_result_file} | sed '1d' | awk '$NF!=\"NOT\"{{ print $1 }}' > {de_mirna_list}",
            shell=True, check=True
        )
        cmd = f"Rscript {dotplot_script} -i {tmp_exp} -o {output_dir} -f {de_mirna_list}"
    elif args.miRNA:
        cmd = f"Rscript {dotplot_script} -i {tmp_exp} -o {output_dir} --miRNA \"{args.miRNA}\""
    elif args.file:
        mir_file = Path(args.file).resolve()
        if not mir_file.is_file():
            sys.exit(f"miRNA list file not found: {mir_file}")
        cmd = f"Rscript {dotplot_script} -i {tmp_exp} -o {output_dir} -f {mir_file}"
    else:
        cmd = f"Rscript {dotplot_script} -i {tmp_exp} -o {output_dir}"

    subprocess.run(cmd, shell=True, check=True)
    print("Expression dotplot generated.")

    # ---- 9. miRNA expression ----
    DEG_file = next(output_dir.glob("*_res.txt"), None)
    if DEG_file is None:
        raise FileNotFoundError(f"No *_res.txt file found in {output_dir}")
    print(f"Found DEG file: {DEG_file}")

    miRNAexp_script = scripts_dir / "miRNA_expression.R"
    cmd = (f"Rscript {miRNAexp_script} -i {tmp_exp} --deg {DEG_file} "
           f"-o {output_dir}")
    subprocess.run(cmd, shell=True, check=True)
    print("miRNA expression analysis completed.")

    # ---- 10. Create final expression file ----
    DEG_file = next(output_dir.glob("*_res.txt"), None)
    if DEG_file is None:
        raise FileNotFoundError(f"No *_res.txt file found in {output_dir}")
    print(f"Found DEG file: {DEG_file}")

    cmd_header = (
        f"printf 'miRNA\\t' > {output_dir}/final_miRNA_expression.txt; "
        f"head -1 {DEG_file} | cut -f2- | tr '\\n' '\\t' >> {output_dir}/final_miRNA_expression.txt; "
        f"head -1 {tmp_exp} | cut -f2- >> {output_dir}/final_miRNA_expression.txt"
    )
    subprocess.run(cmd_header, shell=True, check=True, executable='/bin/bash')

    cmd_join = (
        f"join -t $'\\t' -1 1 -2 1 "
        f"<(tail -n +2 {DEG_file} | sort -k1,1) "
        f"<(tail -n +2 {tmp_exp} | sort -k1,1) >> {output_dir}/final_miRNA_expression.txt"
    )
    subprocess.run(cmd_join, shell=True, check=True, executable='/bin/bash')
    print(f"Final file created: {output_dir}/final_miRNA_expression.txt")

    # ---- 11. Cleanup ----
    shutil.rmtree(temp_dir)
    print(f"\nDifferential expression analysis completed. Results in {output_dir}")