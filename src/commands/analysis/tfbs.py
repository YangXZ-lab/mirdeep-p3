#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
TFBS: Transcription factor binding site analysis (PlantTFDB based).
Supports input via -i (basic-info) or -b (BED file).
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
    """Define TFBS-specific arguments."""
    parser.add_argument("-i", "--input",
                        help="Input basic-info file(s) from annotation, comma separated")
    parser.add_argument("-b", "--bed",
                        help="Input BED file (instead of basic-info)")
    parser.add_argument("-o", "--output",
                        help="Output directory (default: mirdeep-tfbs-<timestamp>)")
    parser.add_argument("--fai",
                        help="Genome FASTA index file (.fai) [required]")
    parser.add_argument("--list", action="store_true",
                        help="List available species from PlantTFDB and exit")
    parser.add_argument("-s", "--species", default="Arabidopsis_thaliana",
                        help="Species name, quoted (default: Arabidopsis_thaliana)")
    parser.add_argument("-u", "--upstream", type=int, default=2000,
                        help="Upstream length to extract (default: 2000)")
    parser.add_argument("-e", "--evalue", default="1e-6",
                        help="E-value threshold for FIMO (default: 1e-6)")
    parser.add_argument("-g", "--genome",
                        help="Reference genome FASTA file [required]")
    parser.add_argument("-p", "--picture", action="store_true",
                        help="Generate TFBS report picture")


def run(args):
    """Execute TFBS analysis."""
    # Resolve project root
    project_root = getattr(args, 'project_root', None) or Path(__file__).resolve().parents[2]
    motif_dir = project_root / "data" / "motif"

    # ---- Handle --list ----
    if args.list:
        sp_map_path = motif_dir / "sp.map"
        if not sp_map_path.is_file():
            sys.exit(f"Species map file not found: {sp_map_path}")
        print("Available species for -s/--species:\n")
        with open(sp_map_path) as f:
            for line in f:
                line = line.strip()
                if line:
                    cols = line.split('\t')
                    if cols:
                        print(f"  {cols[0]}")
        sys.exit(0)


    if not args.fai:
        sys.exit("Error: --fai is required. Use -h for help.")
    if not args.genome:
        sys.exit("Error: -g/--genome is required. Use -h for help.")
        
    # ---- Check required external tools ----
    tools_ok, missing = check_external_tools({'bedtools': None, 'fimo': None})
    if not tools_ok:
        sys.exit("Error: missing required external dependencies (bedtools, fimo).")

    # ---- Determine input ----
    if args.input and args.bed:
        sys.exit("Error: please specify either -i (basic-info) or -b (BED), not both.")
    if not args.input and not args.bed:
        sys.exit("Error: you must provide either -i or -b.")

    # ---- Output directory ----
    if args.output:
        output_dir = Path(args.output)
    else:
        output_dir = Path(f"mirdeep-tfbs-{datetime.now().strftime('%m%d%y-%H%M')}")
    output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = output_dir / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)

    input_bed = temp_dir / "tfbs_input.bed"

    # ---- Step 1: Prepare BED file ----
    if args.input:
        basic_files = [Path(f.strip()) for f in args.input.split(",")]
        with open(input_bed, 'w') as outf:
            for bf in basic_files:
                if not bf.is_file():
                    sys.exit(f"Basic-info file not found: {bf}")
                subprocess.run(
                    f"cat {bf} | awk '{{ print $5\"\\t\"$7\"\\t\"$8\"\\t\"$1\"\\t\"$2\"\\t\"$6 }}' >> {input_bed}",
                    shell=True, check=True
                )
    else:
        bed_file = Path(args.bed)
        if not bed_file.is_file():
            sys.exit(f"BED file not found: {bed_file}")
        shutil.copy2(bed_file, input_bed)

    # ---- Step 2: Upstream extraction ----
    upstream_bed = temp_dir / "tfbs_input.upstream.bed"
    get_upstream_script = project_root / "scripts" / "get_upstream.py"
    cmd = (f"python {get_upstream_script} -u {args.upstream} "
           f"-o {upstream_bed} {input_bed} {args.fai}")
    subprocess.run(cmd, shell=True, check=True)

    # ---- Step 3: Get fasta sequences ----
    upstream_fa = temp_dir / "tfbs_input.upstream.fa"
    genome_file = Path(args.genome)
    if not genome_file.is_file():
        sys.exit(f"Genome file not found: {genome_file}")
    subprocess.run(
        f"bedtools getfasta -s -nameOnly -bed {upstream_bed} "
        f"-fi {genome_file} -fo {upstream_fa}",
        shell=True, check=True
    )

    # ---- Step 4: Lookup species motif files ----
    sp_map_file = motif_dir / "sp.map"
    if not sp_map_file.is_file():
        sys.exit(f"Species map file not found: {sp_map_file}")
    species = args.species.strip()
    motif_files = None
    with open(sp_map_file) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                sys.exit(f"Error: sp.map line {line_num} has insufficient columns (expected 4).")
            if parts[0].strip() == species:
                motif_files = {
                    'tf_list': parts[1].strip(),
                    'meme': parts[2].strip(),
                    'bg': parts[3].strip()
                }
                break
    if motif_files is None:
        sys.exit(f"Species '{species}' not found in {sp_map_file}. Use --list to see available options.")

    tf_list_file = motif_dir / motif_files['tf_list']
    meme_file = motif_dir / motif_files['meme']
    bg_file = motif_dir / motif_files['bg']
    for fpath in (tf_list_file, meme_file, bg_file):
        if not fpath.is_file():
            sys.exit(f"Motif file not found: {fpath}")

    # ---- Step 5: Run FIMO ----
    fimo_out = temp_dir
    log_file = output_dir / "tfbs.log"
    with open(log_file, 'w') as log_fh:
        subprocess.run(
            f"fimo --max-strand --bgfile {bg_file} --thresh {args.evalue} "
            f"-oc {fimo_out} {meme_file} {upstream_fa}",
            shell=True, check=True, stdout=log_fh, stderr=subprocess.STDOUT
        )

    # ---- Step 6: Post‑process FIMO output ----
    tfbs_output = output_dir / "tfbs.tsv"
    subprocess.run(
        f"awk 'NR==FNR {{a[$2]=$3;next}} {{if($1 in a){{ print $0\"\\t\"a[$1] }}}}' "
        f"{tf_list_file} {fimo_out}/fimo.tsv | "
        f"awk '{{ print $1\"\\t\"$NF\"\\t\"$3\"\\t\"$4\"-\"$5\"\\t\"$6\"\\t\"$8\"\\t\"$9\"\\t\"$10 }}' | "
        f"awk 'BEGIN{{FS=OFS=\"\\t\"}} {{gsub(/\\([+-]\\)/,\"\",$3); print}}' > {tfbs_output}",
        shell=True, check=True
    )

    # ---- Step 7: Optional picture ----
    if args.picture:
        rscript_script = project_root / "scripts" / "tfbs_report.R"
        if shutil.which("Rscript") is None:
            print("Warning: Rscript not found, skipping picture generation.")
        else:
            subprocess.run(
                f"Rscript {rscript_script} -i {tfbs_output} -o {output_dir}",
                shell=True, check=True
            )

    # Clean temp
    shutil.rmtree(temp_dir)

    print(f"\nTFBS analysis completed. Results written to {tfbs_output}")
