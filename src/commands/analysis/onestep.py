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
from utils.shellquote import shq


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
    parser.add_argument("--EGGNOG_DATA_DIR",
                        help="Path to a local eggNOG database directory. Passed through to "
                             "Functional_analysis as emapper.py --data_dir. Without it, "
                             "eggNOG-mapper falls back to its own default database location "
                             "and may try to download it.")
    parser.add_argument("--kojson",
                        help="Path to ko00001.json (KO hierarchy used to build the OrgDb). "
                             "Default: <project_root>/data/ko00001.json")


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

    # ---- eggNOG database directory ----
    # Checked here, before the output directory is erased and before the four
    # earlier steps run: a bad path would otherwise only surface at the very end,
    # after the whole pipeline had been wasted.
    eggnog_data_dir = getattr(args, "EGGNOG_DATA_DIR", None)
    if eggnog_data_dir:
        eggnog_path = Path(eggnog_data_dir).expanduser()
        if not eggnog_path.is_dir():
            sys.exit(f"Error: EGGNOG_DATA_DIR not found: {eggnog_path}")
        try:
            empty = not any(eggnog_path.iterdir())
        except OSError:
            empty = False
        if empty:
            print(f"[onestep] WARNING: eggNOG data directory is empty: {eggnog_path}")
        else:
            print(f"[onestep] Using eggNOG database: {eggnog_path.resolve()}")
        args.EGGNOG_DATA_DIR = str(eggnog_path)

    # ---- ko00001.json (same reasoning: fail here, not after five wasted steps) ----
    kojson = getattr(args, "kojson", None)
    if kojson:
        kojson_path = Path(kojson).expanduser()
        if not kojson_path.is_file():
            sys.exit(f"Error: ko00001.json not found: {kojson_path}")
        print(f"[onestep] Using KO hierarchy: {kojson_path.resolve()}")
        args.kojson = str(kojson_path)
    
    output_dir = Path(args.output)
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- Job-local R library ----
    # The OrgDb is rebuilt on every run, so it must not be installed into the
    # shared conda R library: concurrent jobs would deadlock on 00LOCK-* and the
    # base environment would be polluted.  Point every R subprocess at a per-job
    # library instead; the variable is inherited by all child processes.
    r_lib = output_dir / "Rlib"
    r_lib.mkdir(parents=True, exist_ok=True)
    os.environ["MIRDEEP_R_LIB"] = str(r_lib)
    print(f"[onestep] Using job-local R library: {r_lib}")

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
            f"cat {shq(deg_res_file)} | sed '1d' | awk '$NF!=\"NOT\"{{ print $1 }}' | "
            f"awk 'NR==FNR {{a[$1];next}} {{if($1 in a){{ print $0 }}}}' - "
            f"{shq(output_dir / 'target_finder.tsv')} > {shq(target_de)}"
        )
        subprocess.run(awk_cmd, shell=True, check=True)
        target_input = str(target_de)
    else:
        target_input = str(output_dir / "target_finder.tsv")

    functional_partial = False
    func_cmd = [python, str(main_script), "analysis", "Functional_analysis",
                "-p", args.protein,
                "--target", target_input,
                "-t", str(args.threads),
                "-o", str(output_dir)]
    if args.chord:
        func_cmd.append("--chord")
    if getattr(args, "EGGNOG_DATA_DIR", None):
        func_cmd += ["--EGGNOG_DATA_DIR", str(args.EGGNOG_DATA_DIR)]
    if getattr(args, "kojson", None):
        func_cmd += ["--kojson", str(args.kojson)]
    print("[onestep] Running Functional_analysis...")
    # 0 ok | 2 PARTIAL_SUCCESS | 3 NO_TERMS (empty result) | other = failure.
    # check=True would turn 2 and 3 into a traceback; classify them instead.
    func_rc = subprocess.run(func_cmd).returncode
    if func_rc == 2:
        functional_partial = True
        print("[onestep] WARNING: Functional_analysis reported PARTIAL_SUCCESS -- "
              "some enrichment sub-tasks failed, see its summary above.")
    elif func_rc == 3:
        print("[onestep] NOTE: Functional_analysis found no significant enrichment "
              "terms (empty result, not an error).")
    elif func_rc != 0:
        raise RuntimeError(f"Functional_analysis failed (exit {func_rc}); "
                           f"see its output above")

    # ---- Cleanup temp ----
    temp_dir = output_dir / "temp"
    if temp_dir.exists():
        shutil.rmtree(temp_dir)

    if functional_partial:
        print("Onestep pipeline finished with PARTIAL_SUCCESS: every step ran, but some "
              "enrichment sub-tasks failed -- this run is NOT a clean success and its "
              "enrichment output is incomplete.")
        sys.exit(2)
    print("Onestep pipeline completed successfully.")