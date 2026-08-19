#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
PmiREN annotation pipeline: BLAST‑based family assignment and renaming.

This script automates the following steps:
  1. Build BLAST database for the reference miRNA family (optional, or use provided index)
  2. Align query sequences against the reference; keep alignments >= 13 bp
  3. Score, filter, and assign families to aligned sequences
  4. Extract unaligned sequences and run a self‑alignment/clustering step
  5. Produce a unified family mapping and rename the query FASTA accordingly

Output:
  - output/anno.fasta     : renamed sequences
  - output/anno.map       : updated mapping file (original → new)
  - output/temp/          : all intermediate files
  - output/index/         : BLAST databases (unless -d is used for the first DB)

Usage:
  python build_annotation.py -i query.fasta -p reference.fasta -o outdir [options]
"""

import argparse
import os
import sys
import subprocess
import shutil
from pathlib import Path


def run_cmd(cmd, description="", exit_on_fail=True):
    """Print and run a shell command. Exit if it fails and exit_on_fail is True."""
    if description:
        sys.stderr.write(f"[STEP] {description}\n")
    sys.stderr.write(f"[CMD] {cmd}\n")
    result = subprocess.run(cmd, shell=True, executable='/bin/bash')
    if exit_on_fail and result.returncode != 0:
        sys.stderr.write(f"Error: command failed with exit code {result.returncode}\n")
        sys.exit(1)
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Annotate novel miRNAs and rename them according to a reference database.")
    parser.add_argument('-i', '--input', required=True,
                        help="Input FASTA file with query sequences")
    parser.add_argument('-p', '--pmiren', required=True,
                        help="Reference FASTA (core dataset)")
    parser.add_argument('-o', '--output', required=True,
                        help="Output directory")
    parser.add_argument('-d', '--index', default=None,
                        help="Pre-built BLAST database prefix (skip building the first DB)")
    parser.add_argument('--threads', type=int, default=1,
                        help="Number of threads for BLAST (default: 1)")
    parser.add_argument('-t', '--threshold', type=float, default=70.0,
                        help="Score threshold for filtering (default: 70)")
    parser.add_argument('--type', choices=['MIR', 'MIRN'], default='MIR',
                        help="Prefix type for newly assigned families (default: MIR)")
    parser.add_argument('-s', '--start', type=int, default=None,
                        help="Starting family number (required when --type MIR)")
    parser.add_argument('--prefix', default=None,
                        help="Optional prefix for renamed sequences (e.g., 'Ath')")
    args = parser.parse_args()

    # ---- Determine project root and scripts directory ----
    project_root = Path(__file__).resolve().parent.parent
    scripts_dir = project_root / 'scripts'
    sys.stderr.write(f"Project root: {project_root}\n")

    # ---- Validate arguments ----
    if args.type == 'MIR' and args.start is None:
        sys.exit("Error: -s/--start is required when --type MIR is used.")
    if args.type == 'MIRN' and args.start is not None:
        sys.exit("Error: -s/--start cannot be used with --type MIRN.")

    input_fa = os.path.abspath(args.input)
    pmiren_fa = os.path.abspath(args.pmiren)
    output_dir = os.path.abspath(args.output)

    # ---- Create output directories ----
    temp_dir = os.path.join(output_dir, 'temp')
    index_dir = os.path.join(output_dir, 'index')
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(index_dir, exist_ok=True)

    # ---- Build or use reference BLAST database ----
    if args.index:
        ref_db = os.path.abspath(args.index)
        sys.stderr.write(f"Using provided reference DB: {ref_db}\n")
    else:
        ref_db = os.path.join(index_dir, 'pmiren_index')
        run_cmd(f"makeblastdb -in {pmiren_fa} -dbtype nucl -out {ref_db}",
                "Building BLAST database for reference sequences")

    # ---- Define file paths for intermediate outputs ----
    anno_aln = os.path.join(temp_dir, 'anno.aln')
    anno_aln_filter = os.path.join(temp_dir, 'anno.aln.filter')
    anno_total = os.path.join(temp_dir, 'anno_total.fasta')
    anno_score = os.path.join(temp_dir, 'anno.aln.filter.score')
    anno_score_stat = os.path.join(temp_dir, 'anno.aln.filter.score.stat')
    anno_score_filter = os.path.join(temp_dir, 'anno.aln.filter.score.stat.filter')
    anno_assigned = os.path.join(temp_dir, 'anno.aln.filter.score.stat.filter.anno')
    anno_assigned_stat = os.path.join(temp_dir, 'anno.aln.filter.score.stat.filter.stat')

    nonaln_fa = os.path.join(temp_dir, 'anno_nonaln.fasta')
    nonaln_db = os.path.join(index_dir, 'anno_nonaln')
    nonaln_aln = os.path.join(temp_dir, 'anno_nonaln.aln')
    nonaln_aln_filter = os.path.join(temp_dir, 'anno_nonaln.aln.filter')
    nonaln_score = os.path.join(temp_dir, 'anno_nonaln.aln.filter.score')
    single_file = os.path.join(temp_dir, 'anno_nonaln.aln.filter.score.single')
    multi_file = os.path.join(temp_dir, 'anno_nonaln.aln.filter.score.multi')
    multi_score = os.path.join(temp_dir, 'anno_nonaln.aln.filter.score.multi.score')
    multi_score_filter = os.path.join(temp_dir, 'anno_nonaln.aln.filter.score.multi.score.filter')
    cluster_file = os.path.join(temp_dir, 'anno_nonaln.aln.filter.score.multi.score.cl')
    anno_map = os.path.join(temp_dir, 'anno.map')

    final_fasta = os.path.join(output_dir, 'anno.fasta')
    final_map = os.path.join(output_dir, 'anno.map')

    # ---- External tool/script paths ----
    scoring_script = scripts_dir / 'bowtie_scoring_blastn.py'
    stat_script = scripts_dir / 'blast-score-stat.py'
    filter_script = scripts_dir / 'blast-score-filter.py'
    assign_script = scripts_dir / 'assign_family_by_score.py'
    bowtie_filter_script = scripts_dir / 'bowtie_score_filter.py'
    cluster_script = scripts_dir / 'MIR_cluster.py'
    remap_script = scripts_dir / 'reanno_remap_v5.py'
    rename_script = scripts_dir / 'reanno_fasta_by_map_v3.py'

    # ================================================================
    # STEP 1: BLAST alignment of query against reference
    # ================================================================
    run_cmd(
        f"blastn -task blastn-short "
        f"-query {input_fa} "
        f"-db {ref_db} "
        f"-outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue qseq sseq\" "
        f"-out {anno_aln} "
        f"-num_threads {args.threads}",
        "Aligning query sequences against reference database"
    )

    run_cmd(f"cat {anno_aln} | awk '$4>=13{{ print $0 }}' > {anno_aln_filter}",
            "Filtering alignments with length >= 13")

    run_cmd(f"cat {input_fa} {pmiren_fa} > {anno_total}",
            "Merging input and reference FASTA for scoring")

    run_cmd(f"python {scoring_script} -i {anno_aln_filter} -f {anno_total} -o {anno_score}",
            "Scoring filtered alignments")

    run_cmd(f"python {stat_script} -i {anno_score} -o {anno_score_stat}",
            "Extracting per-query best scores")

    run_cmd(f"python {filter_script} -s {anno_score} -m {anno_score_stat} -t {args.threshold} -o {anno_score_filter}",
            f"Applying score threshold (>= {args.threshold})")

    run_cmd(
        f"python {assign_script} -i {anno_score_filter} "
        f"-f {input_fa} "
        f"-o {anno_assigned} "
        f"-s {anno_assigned_stat}",
        "Assigning families to aligned queries"
    )

    # ================================================================
    # STEP 2: Extract unaligned sequences
    # ================================================================
    run_cmd(
        f"cat {input_fa} | grep '>' | sed 's/>//g' | "
        f"awk 'NR==FNR {{a[$1];next}} {{if(!($1 in a)){{ print $0 }}}}' {anno_assigned} - | "
        f"seqkit grep -n -f - {input_fa} -o {nonaln_fa}",
        "Extracting unassigned sequences"
    )

    # ================================================================
    # STEP 3: Self-alignment and clustering of unaligned sequences
    # ================================================================
    run_cmd(f"makeblastdb -in {nonaln_fa} -dbtype nucl -out {nonaln_db}",
            "Building BLAST database for unassigned sequences")

    run_cmd(
        f"blastn -task blastn-short "
        f"-query {nonaln_fa} "
        f"-db {nonaln_db} "
        f"-outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue qseq sseq\" "
        f"-out {nonaln_aln} "
        f"-num_threads {args.threads}",
        "Self-aligning unassigned sequences"
    )

    run_cmd(f"cat {nonaln_aln} | awk '$4>=13{{ print $0 }}' > {nonaln_aln_filter}",
            "Filtering self-alignments with length >= 13")

    run_cmd(f"python {scoring_script} -i {nonaln_aln_filter} -f {nonaln_fa} -o {nonaln_score}",
            "Scoring self-alignments")

    # Single and multi classification
    run_cmd(
        f"cat {nonaln_score} | sed '1d' | awk '$NF>=75{{ print $0 }}' | "
        f"awk '{{sum[$1]+=1;a[$2]+=1}}END{{for(i in sum){{ print i\"\\t\"sum[i]\"\\t\"a[i] }}}}' | "
        f"awk '$2==1&&$3==1{{ print $1 }}' > {single_file}",
        "Identifying single-copy unassigned miRNAs"
    )

    run_cmd(
        f"cat {nonaln_score} | sed '1d' | awk '$NF>=75{{ print $0 }}' | "
        f"awk '{{sum[$1]+=1;a[$2]+=1}}END{{for(i in sum){{ print i\"\\t\"sum[i]\"\\t\"a[i] }}}}' | "
        f"awk '$2!=1||$3!=1{{ print $1 }}' > {multi_file}",
        "Identifying multi-copy unassigned miRNAs"
    )

    # Extract multi scores (non-self, both query and target in multi set)
    run_cmd(
        f"cat {nonaln_score} | "
        f"awk 'NR==FNR {{a[$1];next}} {{if($1 in a){{ print $0 }}}}' {multi_file} - | "
        f"awk '$1!=$2{{ print $0 }}' | "
        f"awk 'NR==FNR {{a[$1];next}} {{if($2 in a){{ print $0 }}}}' {multi_file} - > {multi_score}",
        "Extracting multi-copy alignments for clustering"
    )

    run_cmd(f"python {bowtie_filter_script} -i {multi_score} -o {multi_score_filter}",
            "Filtering multi-copy scores")

    run_cmd(f"python {cluster_script} -i {multi_score_filter} -o {cluster_file}",
            "Clustering multi-copy miRNAs")

    # ================================================================
    # STEP 4: Merge all families into a single mapping
    # ================================================================
    remap_cmd = (
        f"python {remap_script} "
        f"-f {pmiren_fa} "
        f"--alncl {anno_assigned} "
        f"--single {single_file} "
        f"--multi {cluster_file} "
        f"--type {args.type} "
    )
    if args.type == 'MIR':
        remap_cmd += f"-s {args.start} "
    remap_cmd += f"-o {anno_map}"
    run_cmd(remap_cmd, "Creating unified family mapping")

    # ================================================================
    # STEP 5: Rename query FASTA
    # ================================================================
    rename_cmd = f"python {rename_script} -m {anno_map} -i {input_fa} -o {final_fasta} -r {final_map}"
    if args.prefix:
        rename_cmd += f" --prefix {args.prefix}"
    run_cmd(rename_cmd, "Renaming query sequences according to mapping")

    sys.stderr.write("Pipeline completed successfully.\n")
    sys.stderr.write(f"Renamed FASTA: {final_fasta}\n")
    sys.stderr.write(f"Updated mapping: {final_map}\n")


if __name__ == "__main__":
    main()
