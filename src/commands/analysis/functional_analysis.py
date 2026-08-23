#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
Functional_analysis: miRNA target functional enrichment.
Supports building OrgDb from eggNOG or using pre-built OrgDb,
and handles gene or miRNA-target input.
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
    """Define Functional_analysis-specific arguments."""
    # Input modes
    parser.add_argument("-p", "--protein",
                        help="Protein FASTA file for eggNOG annotation (starts full pipeline)")
    parser.add_argument("--orgdb",
                        help="Path to OrgDb directory (org.Morg.eg.db) if already built")
    parser.add_argument("-f", "--file",
                        help="Path to directory containing pathway2gene and pathway2name files")
    parser.add_argument("-g", "--gene",
                        help="Gene list file for direct enrichment analysis")
    parser.add_argument("--target",
                        help="miRNA-target file (at least two columns: miRNA, target)")

    # eggNOG specific
    parser.add_argument("--EGGNOG_DATA_DIR",
                        help="Path to eggNOG data directory (default: data/eggnog_data_dir)")
    parser.add_argument("--kojson",
                        help="Path to ko00001.json (default: data/ko00001.json)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads for eggNOG mapper (default: 1)")

    # Output and other
    parser.add_argument("-o", "--output",
                        help="Output directory (default: mirdeep-functional_analysis-<timestamp>)")
    parser.add_argument("--chord", action="store_true",
                        help="Generate chord diagram (only with --target)")


def run(args):
    """Execute Functional_analysis."""
    project_root = getattr(args, 'project_root', None) or Path(__file__).resolve().parents[2]
    data_dir = project_root / "data"
    scripts_dir = project_root / "scripts"

    # ---- 1. Parameter validation ----
    build_orgdb = args.protein is not None
    if build_orgdb and (args.orgdb is not None or args.file is not None):
        sys.exit("Error: -p/--protein cannot be combined with --orgdb or -f. Use -p to build from scratch.")
    if not build_orgdb:
        if not args.orgdb or not args.file:
            sys.exit("Error: when not building (no -p), both --orgdb and -f are required.")
    if args.gene and args.target:
        sys.exit("Error: please specify either -g/--gene or --target, not both.")
    if not args.gene and not args.target:
        sys.exit("Error: you must provide an input for enrichment: -g/--gene or --target.")
    if args.chord and not args.target:
        sys.exit("Error: --chord is only valid with --target.")

    # ---- 2. Output directory ----
    if args.output:
        output_dir = Path(args.output)
    else:
        output_dir = Path(f"mirdeep-functional_analysis-{datetime.now().strftime('%m%d%y-%H%M')}")
    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- 3. Build orgdb if requested ----
    if build_orgdb:
        # Check external tools: emapper.py, Rscript
        tools_ok, missing = check_external_tools({'emapper.py': None, 'Rscript': None})
        if not tools_ok:
            sys.exit("Error: missing required external dependencies (emapper.py, Rscript).")

        protein_fasta = Path(args.protein)
        if not protein_fasta.is_file():
            sys.exit(f"Protein FASTA file not found: {protein_fasta}")

        # Resolve kojson (used later in build_orgdb.R)
        kojson_file = Path(args.kojson) if args.kojson else data_dir / "ko00001.json"
        if not kojson_file.is_file():
            sys.exit(f"ko00001.json not found: {kojson_file}")

        # ---- 3a. Run eggNOG mapper ----
        print("Running eggNOG-mapper...")
        emapper_output_base = output_dir

        # Determine whether to pass --data_dir based on user specification
        if args.EGGNOG_DATA_DIR:
            eggnog_data = Path(args.EGGNOG_DATA_DIR)
            if not eggnog_data.is_dir():
                sys.exit(f"EGGNOG_DATA_DIR not found: {eggnog_data}")
            cmd = (f"emapper.py --data_dir {eggnog_data} --cpu {args.threads} "
                   f"-m diamond --override --dbmem "
                   f"-d euk --tax_scope Viridiplantae -i {protein_fasta} "
                   f"-o {emapper_output_base}")
        else:
            # No custom data dir: rely on emapper's default location
            cmd = (f"emapper.py --cpu {args.threads} -m diamond --override --dbmem "
                   f"-d euk --tax_scope Viridiplantae -i {protein_fasta} "
                   f"-o {emapper_output_base}")
        subprocess.run(cmd, shell=True, check=True)

        # ---- 3b. Process eggNOG outputs ----
        eggnog_annot = emapper_output_base.with_name(emapper_output_base.name + ".emapper.annotations")
        go_annot = emapper_output_base / "Go.eggnog.emapper.annotations"
        subprocess.run(
            f"sed '/^##/d' {eggnog_annot} | sed 's/#//g' | "
            f"awk -vFS='\\t' -vOFS='\\t' '{{print $1,$9,$10,$12}}' > {go_annot}",
            shell=True, check=True
        )

        # ---- 3c. Build OrgDb and pathway files ----
        build_script = scripts_dir / "build_orgdb.R"
        orgdb_outdir = emapper_output_base
        cmd = (f"Rscript {build_script} -i {go_annot} --kojson {kojson_file} -o {orgdb_outdir}")
        subprocess.run(cmd, shell=True, check=True)

        # After building, set orgdb and file_paths for later use
        orgdb_path = orgdb_outdir / "org.Morg.eg.db"
        pathway_dir = orgdb_outdir
    else:
        # Use provided paths
        orgdb_path = Path(args.orgdb)
        pathway_dir = Path(args.file)
        if not orgdb_path.is_dir():
            sys.exit(f"OrgDb directory not found: {orgdb_path}")
        if not pathway_dir.is_dir():
            sys.exit(f"Pathway file directory not found: {pathway_dir}")

    # ---- 4. Functional enrichment ----
    if args.gene:
        enrich_script = scripts_dir / "enrich_analysis.R"
        gene_file = Path(args.gene)
        if not gene_file.is_file():
            sys.exit(f"Gene file not found: {gene_file}")
        subprocess.run(
            f"Rscript {enrich_script} -i {orgdb_path} -f {pathway_dir} "
            f"-g {gene_file} --goterm 10 -o {output_dir}",
            shell=True, check=True
        )
    else:
        target_file = Path(args.target)
        if not target_file.is_file():
            sys.exit(f"Target file not found: {target_file}")
        enrich_script = scripts_dir / "miRNA_enrich_analysis.py"
        cmd = (f"python {enrich_script} -i {target_file} --orgdb {orgdb_path} "
               f"-f {pathway_dir} -o {output_dir}")
        subprocess.run(cmd, shell=True, check=True)

        if args.chord:
            chord_script = scripts_dir / "miRNA_chord_type2.R"
            subprocess.run(
                f"Rscript {chord_script} -i {target_file} --orgdb {orgdb_path} "
                f"-f {pathway_dir} -o {output_dir}",
                shell=True, check=True
            )

    print(f"\nFunctional analysis completed. Results in {output_dir}")
