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
from utils.shellquote import shq


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


def run_enrichment_step(cmd: str, label: str) -> str:
    """
    Run one enrichment command and classify its outcome.

    Returns 'ok', 'partial' or 'no_terms'; raises RuntimeError on a real failure.
    The status codes are the ones documented in scripts/miRNA_enrich_analysis.py:
    0 ok, 2 PARTIAL_SUCCESS, 3 NO_TERMS (empty result), anything else = failure.
    """
    rc = subprocess.run(cmd, shell=True).returncode
    if rc == 0:
        return "ok"
    if rc == 2:
        print(f"[warn] {label}: PARTIAL_SUCCESS -- some sub-tasks failed; "
              f"the enrichment output is incomplete (see the summary above).")
        return "partial"
    if rc == 3:
        print(f"[note] {label}: no significant enrichment terms were found "
              f"(empty result, not an error).")
        return "no_terms"
    raise RuntimeError(f"{label} failed (exit {rc}); see the output above")


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

    # Clear existing output directory to avoid stale files
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Job-local R library for the freshly built OrgDb (mirrors onestep.py so a
    # standalone Functional_analysis run behaves the same way).  setdefault keeps
    # whatever the parent driver already exported.
    r_lib = output_dir / "Rlib"
    r_lib.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MIRDEEP_R_LIB", str(r_lib))
    print(f"[functional_analysis] Using R library: {os.environ['MIRDEEP_R_LIB']}")

    # ---- 3. Build orgdb if requested ----
    if build_orgdb:
        # Check external tools: emapper.py, Rscript
        tools_ok, missing = check_external_tools({'emapper.py': None, 'Rscript': None})
        if not tools_ok:
            sys.exit("Error: missing required external dependencies (emapper.py, Rscript).")

        protein_fasta = Path(args.protein)
        if not protein_fasta.is_file():
            sys.exit(f"Protein FASTA file not found: {protein_fasta}")

        # ---- 3a. Preprocess protein FASTA: remove periods from sequences ----
        protein_clean = output_dir / f"{protein_fasta.stem}.clean{protein_fasta.suffix}"
        print(f"Cleaning protein sequences (removing '.') to {protein_clean} ...")
        cmd = f"sed '/^>/! s/\\.//g' {shq(protein_fasta)} > {shq(protein_clean)}"
        subprocess.run(cmd, shell=True, check=True)
        if not protein_clean.is_file() or protein_clean.stat().st_size == 0:
            sys.exit(f"Error: cleaned protein file was not created or is empty: {protein_clean}")

        # Resolve kojson (used later in build_orgdb.R)
        kojson_file = Path(args.kojson) if args.kojson else data_dir / "ko00001.json"
        if not kojson_file.is_file():
            sys.exit(f"ko00001.json not found: {kojson_file}")

        # ---- 3b. Run eggNOG mapper ----
        print("Running eggNOG-mapper...")
        emapper_output_base = output_dir

        # Determine whether to pass --data_dir based on user specification
        if args.EGGNOG_DATA_DIR:
            eggnog_data = Path(args.EGGNOG_DATA_DIR)
            if not eggnog_data.is_dir():
                sys.exit(f"EGGNOG_DATA_DIR not found: {eggnog_data}")
            cmd = (f"emapper.py --data_dir {shq(eggnog_data)} --cpu {args.threads} "
                   f"-m diamond --override --dbmem "
                   f"-d euk --tax_scope Viridiplantae -i {shq(protein_clean)} "
                   f"-o {shq(emapper_output_base)}")
        else:
            cmd = (f"emapper.py --cpu {args.threads} -m diamond --override --dbmem "
                   f"-d euk --tax_scope Viridiplantae -i {shq(protein_clean)} "
                   f"-o {shq(emapper_output_base)}")
        subprocess.run(cmd, shell=True, check=True)

        # ---- 3c. Process eggNOG outputs ----
        eggnog_annot = emapper_output_base.with_name(emapper_output_base.name + ".emapper.annotations")
        go_annot = emapper_output_base / "Go.eggnog.emapper.annotations"
        subprocess.run(
            f"sed '/^##/d' {shq(eggnog_annot)} | sed 's/#//g' | "
            f"awk -vFS='\\t' -vOFS='\\t' '{{print $1,$9,$10,$12}}' > {shq(go_annot)}",
            shell=True, check=True
        )

        # ---- 3d. Build OrgDb and pathway files ----
        build_script = scripts_dir / "build_orgdb.R"
        orgdb_outdir = emapper_output_base
        cmd = (f"Rscript {shq(build_script)} -i {shq(go_annot)} "
               f"--kojson {shq(kojson_file)} -o {shq(orgdb_outdir)}")
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
        # Gene-based enrichment
        enrich_script = scripts_dir / "enrich_analysis.R"
        gene_file = Path(args.gene)
        if not gene_file.is_file():
            sys.exit(f"Gene file not found: {gene_file}")
        step_status = run_enrichment_step(
            f"Rscript {shq(enrich_script)} -i {shq(orgdb_path)} -f {shq(pathway_dir)} "
            f"-g {shq(gene_file)} --goterm 10 -o {shq(output_dir)}",
            "GO/KEGG enrichment (gene list)")
    else:
        # miRNA-target based enrichment
        target_file = Path(args.target)
        if not target_file.is_file():
            sys.exit(f"Target file not found: {target_file}")
        enrich_script = scripts_dir / "miRNA_enrich_analysis.py"
        cmd = (f"python {shq(enrich_script)} -i {shq(target_file)} "
               f"--orgdb {shq(orgdb_path)} -f {shq(pathway_dir)} -o {shq(output_dir)}")
        step_status = run_enrichment_step(cmd, "miRNA target enrichment")

        # Optional chord diagram
        if args.chord:
            chord_script = scripts_dir / "miRNA_chord_type2.R"
            subprocess.run(
                f"Rscript {shq(chord_script)} -i {shq(target_file)} "
                f"--orgdb {shq(orgdb_path)} -f {shq(pathway_dir)} -o {shq(output_dir)}",
                shell=True, check=True
            )

    if step_status == "partial":
        print(f"\nFunctional analysis finished with PARTIAL_SUCCESS: every step ran, "
              f"but some enrichment sub-tasks failed, so the enrichment output is "
              f"INCOMPLETE. Results in {output_dir}")
        sys.exit(2)

    print(f"\nFunctional analysis completed. Results in {output_dir}")
