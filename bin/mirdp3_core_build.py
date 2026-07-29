#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PmiREN-core V2 construction pipeline.

This script automates the workflow of:
  1. Aligning novel isoform sequences (isoform-out-total) against the
     reference isoform database (isoform-in) using BLASTn-short.
  2. Scoring, filtering, and assigning families to the aligned sequences.
  3. Iteratively clustering unaligned sequences via self-alignment and
     graph-based clustering.
  4. Renaming all novel miRNAs according to a unified naming scheme.
  5. Splitting each family into host and guest sets, and finally merging
     the host set with the reference database to produce the updated core.

Usage:
    python build_pmiren_core_v2.py -i isoform-out-total_unique.fa -s 12178 -o output_dir
"""

import argparse
import os
import shutil
import subprocess
import sys
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
        description="Build PmiREN-core V2 by annotating novel miRNA families.")
    parser.add_argument('-p', '--pmiren', default=None,
                        help="Input isoform-in.fa file (default: <project>/data/isoform-in.fa)")
    parser.add_argument('-i', '--input', required=True,
                        help="Input isoform-out-total_unique.fa file (the query sequences)")
    parser.add_argument('-s', '--start', type=int, required=True,
                        help="Starting family number for new miRNA families (MIRxxxx)")
    parser.add_argument('-o', '--output', required=True,
                        help="Output directory (will be created; temp/ and index/ inside it)")
    parser.add_argument('-c', '--clean', action='store_true',
                        help="Clean (remove) output directory if it already exists before running")
    parser.add_argument('--threshold', type=float, default=75.0,
                        help="Score threshold for single/multi separation (default: 75)")
    args = parser.parse_args()

    # ---- Determine project root ----
    project_root = Path(__file__).resolve().parent.parent
    sys.stderr.write(f"Project root: {project_root}\n")
    scripts_dir = project_root / 'scripts'
    data_dir = project_root / 'data'

    # ---- Set paths ----
    pmiren_fa = args.pmiren if args.pmiren else str(data_dir / 'isoform-in.fa')
    input_fa = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.output)

    # ---- Clean if requested ----
    if args.clean and os.path.exists(output_dir):
        shutil.rmtree(output_dir)
        sys.stderr.write(f"Cleaned output directory: {output_dir}\n")

    # ---- Create directories ----
    temp_dir = os.path.join(output_dir, 'temp')
    index_dir = os.path.join(output_dir, 'index')
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(index_dir, exist_ok=True)

    # ---- Define common file paths ----
    # Output files (placed in output_dir root)
    host_fa = os.path.join(output_dir, 'isoform-out-total-rename-host.fa')
    guest_fa = os.path.join(output_dir, 'isoform-out-total-rename-guest.fa')
    final_in_v2 = os.path.join(output_dir, 'isoform-in-v2.fa')

    # Temporary files
    merged_fa = os.path.join(temp_dir, 'isoform-v2-total.fa')        # combined in+out
    out_aln = os.path.join(temp_dir, 'isoform-out-total.aln')
    out_aln_filter = os.path.join(temp_dir, 'isoform-out-total.aln.filter')
    out_aln_score = os.path.join(temp_dir, 'isoform-out-total.aln.filter.score')
    out_aln_score_stat = os.path.join(temp_dir, 'isoform-out-total.aln.filter.score.stat')
    out_aln_score_stat_filter = os.path.join(temp_dir, 'isoform-out-total.aln.filter.score.stat.filter')
    out_aln_anno = os.path.join(temp_dir, 'isoform-out-total.aln.filter.score.stat.filter.anno')
    out_aln_anno_stat = os.path.join(temp_dir, 'isoform-out-total.aln.filter.score.stat.filter.stat')

    nonaln_fa = os.path.join(temp_dir, 'isoform-out-total-nonaln.fa')
    nonaln_aln = os.path.join(temp_dir, 'isoform-out-total-nonaln.aln')
    nonaln_aln_filter = os.path.join(temp_dir, 'isoform-out-total-nonaln.aln.filter')
    nonaln_aln_score = os.path.join(temp_dir, 'isoform-out-total-nonaln.aln.filter.score')
    nonaln_single = os.path.join(temp_dir, 'isoform-out-total-nonaln.aln.filter.score.single')
    nonaln_multi = os.path.join(temp_dir, 'isoform-out-total-nonaln.aln.filter.score.multi')
    nonaln_multi_score = os.path.join(temp_dir, 'isoform-out-total-nonaln.aln.filter.score.multi.score')
    nonaln_multi_score_filter = os.path.join(temp_dir, 'isoform-out-total-nonaln.aln.filter.score.multi.score.filter')
    nonaln_multi_cl = os.path.join(temp_dir, 'isoform-out-total-nonaln.aln.filter.score.multi.score.cl')

    rename_map = os.path.join(temp_dir, 'isoform-out-total.map')
    rename_fa = os.path.join(temp_dir, 'isoform-out-total-rename.fa')
    rename_updated_map = os.path.join(temp_dir, 'isoform-out-total-rename.map')

    # ---- External tools ----
    # Assumes blastn, makeblastdb, seqkit, awk etc. are in PATH
    scoring_script = str(scripts_dir / 'bowtie_scoring_blastn.py')
    # Scripts in scripts_dir (as specified by user)
    stat_script = str(scripts_dir / 'blast-score-stat.py')
    filter_script = str(scripts_dir / 'blast-score-filter.py')
    filter_script_2 = str(scripts_dir / 'bowtie_score_filter.py')
    assign_script = str(scripts_dir / 'assign_family_by_score.py')
    cluster_script = str(scripts_dir / 'MIR_cluster.py')
    remap_script = str(scripts_dir / 'reanno_remap_v3.py')
    rename_fa_script = str(scripts_dir / 'reanno_fasta_by_map_v2.py')
    split_script = str(scripts_dir / 'split_fasta_gh_v2.py')

    # ===================================================================
    # STEP 1: BLAST alignment of out vs in
    # ===================================================================
    run_cmd(f"makeblastdb -in {pmiren_fa} -dbtype nucl -out {index_dir}/isoform-in",
            "Building BLAST database for isoform-in")

    run_cmd(
        f"blastn -task blastn-short "
        f"-query {input_fa} "
        f"-db {index_dir}/isoform-in "
        f"-outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue qseq sseq\" "
        f"-out {out_aln} "
        f"-num_threads 14",
        "Aligning isoform-out-total against isoform-in (BLASTn-short)"
    )

    # ===================================================================
    # STEP 2: Pre‑filter (length >= 16) and score
    # ===================================================================
    run_cmd(f"cat {out_aln} | awk '$4>=16{{ print $0 }}' > {out_aln_filter}",
            "Filtering alignments with length >= 16")

    run_cmd(f"cat {pmiren_fa} {input_fa} > {merged_fa}",
            "Merging isoform-in and isoform-out-total for scoring")

    run_cmd(f"python {scoring_script} -i {out_aln_filter} -f {merged_fa} -o {out_aln_score}",
            "Scoring filtered alignments")

    # ===================================================================
    # STEP 3: Per‑query maximum and global filter
    # ===================================================================
    run_cmd(f"python {stat_script} -i {out_aln_score} -o {out_aln_score_stat}",
            "Extracting per‑query best score")

    run_cmd(f"python {filter_script} -s {out_aln_score} -m {out_aln_score_stat} -o {out_aln_score_stat_filter}",
            "Applying dynamic score threshold to alignments")

    # ===================================================================
    # STEP 4: Assign family to each matched query
    # ===================================================================
    run_cmd(
        f"python {assign_script} -i {out_aln_score_stat_filter} "
        f"-f {input_fa} "
        f"-o {out_aln_anno} "
        f"-s {out_aln_anno_stat}",
        "Assigning miRNA families to aligned queries"
    )

    # ===================================================================
    # STEP 5: Extract unaligned sequences
    # ===================================================================
    run_cmd(
        f"cat {input_fa} | grep '>' | sed 's/>//g' | "
        f"awk 'NR==FNR {{a[$1];next}} {{if(!($1 in a)){{ print $0 }}}}' "
        f"{out_aln_anno} - | "
        f"seqkit grep -n -f - {input_fa} -o {nonaln_fa}",
        "Extracting sequences that were not assigned to any family"
    )

    # ===================================================================
    # STEP 6: Self‑alignment of unaligned sequences
    # ===================================================================
    run_cmd(f"makeblastdb -in {nonaln_fa} -dbtype nucl -out {index_dir}/isoform-out-nonaln",
            "Building BLAST database for unaligned sequences")

    run_cmd(
        f"blastn -task blastn-short "
        f"-query {nonaln_fa} "
        f"-db {index_dir}/isoform-out-nonaln "
        f"-outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue qseq sseq\" "
        f"-out {nonaln_aln} "
        f"-num_threads 14",
        "Self‑aligning unaligned sequences (BLASTn-short)"
    )

    # ===================================================================
    # STEP 7: Filter, score and separate single vs multi
    # ===================================================================
    run_cmd(f"cat {nonaln_aln} | awk '$4>=16{{ print $0 }}' > {nonaln_aln_filter}",
            "Filtering self‑alignments with length >= 16")

    run_cmd(f"python {scoring_script} -i {nonaln_aln_filter} -f {nonaln_fa} -o {nonaln_aln_score}",
            "Scoring self‑alignments")

    # Single copy miRNAs (score >= threshold, only one hit, one target)
    run_cmd(
        f"cat {nonaln_aln_score} | sed '1d' | awk '$NF>={args.threshold}{{ print $0 }}' | "
        f"awk '{{sum[$1]+=1;a[$2]+=1}}END{{for(i in sum){{ print i\"\\t\"sum[i]\"\\t\"a[i] }}}}' | "
        f"awk '$2==1&&$3==1{{ print $1 }}' > {nonaln_single}",
        f"Identifying single‑copy miRNAs (threshold >= {args.threshold})"
    )

    # Multi‑copy miRNAs (everything else above threshold)
    run_cmd(
        f"cat {nonaln_aln_score} | sed '1d' | awk '$NF>={args.threshold}{{ print $0 }}' | "
        f"awk '{{sum[$1]+=1;a[$2]+=1}}END{{for(i in sum){{ print i\"\\t\"sum[i]\"\\t\"a[i] }}}}' | "
        f"awk '$2!=1||$3!=1{{ print $1 }}' > {nonaln_multi}",
        f"Identifying multi‑copy miRNAs (threshold >= {args.threshold})"
    )

    # ===================================================================
    # STEP 8: Process multi‑copy miRNAs (filter, cluster)
    # ===================================================================
    run_cmd(
        f"cat {nonaln_aln_score} | "
        f"awk 'NR==FNR {{a[$1];next}} {{if($1 in a){{ print $0 }}}}' {nonaln_multi} - | "
        f"awk '$1!=$2{{ print $0 }}' | "
        f"awk 'NR==FNR {{a[$1];next}} {{if($2 in a){{ print $0 }}}}' {nonaln_multi} - > {nonaln_multi_score}",
        "Extracting multi‑copy alignments (non‑self, mutual)"
    )

    run_cmd(f"python {filter_script_2} -i {nonaln_multi_score} -o {nonaln_multi_score_filter}",
            "Filtering multi‑copy scores (keeping best per query‑target pair)")

    run_cmd(f"python {cluster_script} -i {nonaln_multi_score_filter} -o {nonaln_multi_cl}",
            "Clustering multi‑copy miRNAs")

    # ===================================================================
    # STEP 9: Build unified family mapping
    # ===================================================================
    run_cmd(
        f"python {remap_script} --alncl {out_aln_anno} "
        f"--single {nonaln_single} "
        f"--multi {nonaln_multi_cl} "
        f"-s {args.start} "
        f"-o {rename_map}",
        "Merging all family assignments into a single map"
    )

    # ===================================================================
    # STEP 10: Rename sequences in the main input FASTA
    # ===================================================================
    run_cmd(
        f"python {rename_fa_script} -m {rename_map} "
        f"-i {input_fa} "
        f"-o {rename_fa} "
        f"-r {rename_updated_map}",
        "Renaming miRNA sequences according to the unified mapping"
    )

    # ===================================================================
    # STEP 11: Split families into host and guest
    # ===================================================================
    run_cmd(
        f"python {split_script} -m {rename_updated_map} "
        f"-i {rename_fa} "
        f"-s {args.start} "
        f"--host {host_fa} "
        f"--guest {guest_fa}",
        "Splitting renamed families into host and guest sets"
    )

    # ===================================================================
    # STEP 12: Final assembly of the new isoform-in core
    # ===================================================================
    run_cmd(f"cat {pmiren_fa} {host_fa} > {final_in_v2}",
            "Creating the updated isoform-in-v2.fa database")

    sys.stderr.write("Pipeline completed successfully.\n")
    sys.stderr.write(f"Output files:\n  {final_in_v2}\n  {host_fa}\n  {guest_fa}\n")


if __name__ == "__main__":
    main()
