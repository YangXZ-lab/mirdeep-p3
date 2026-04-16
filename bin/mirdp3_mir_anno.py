#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.


"""
mirdp3_mir_anno.py

Main pipeline for miRNA isoform annotation and renaming.

Usage:
    python mirdp3_mir_anno.py -p PmiREN-core.fa -i input.fa -s 11100 [-d index_dir] [-o out_dir] [-t temp_dir] [-clean]
"""

import sys
import os
import re
import argparse
import subprocess
import shutil
import shlex
from datetime import datetime

def check_dependencies():
    """Check if required tools (bowtie, seqkit) are available in PATH."""
    missing = []
    for tool in ["bowtie", "seqkit"]:
        if shutil.which(tool) is None:
            missing.append(tool)
    if missing:
        sys.stderr.write(f"Error: missing dependencies: {', '.join(missing)}\n")
        sys.exit(1)
    else:
        print("[✓] bowtie found")
        print("[✓] seqkit found")

def ensure_dir(path):
    """Create directory if it does not exist."""
    if not os.path.exists(path):
        os.makedirs(path)
    return path

def check_index_exists(index_prefix):
    """Check if Bowtie index files exist (e.g., index_prefix.1.ebwt)."""
    return os.path.exists(f"{index_prefix}.1.ebwt")

def main():
    # Record full command line for logging
    full_command = " ".join(sys.argv)

    parser = argparse.ArgumentParser(description="MiRDP3 mirna annotation pipeline")
    parser.add_argument("-p", "--pmiren", required=True, help="PmiREN-core dataset (FASTA)")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file (isoforms)")
    parser.add_argument("-s", "--start", required=True, type=int, help="Start number for renaming")
    parser.add_argument("-d", "--index", help="Existing Bowtie index prefix (e.g., /path/to/isoform-in)")
    parser.add_argument("-o", "--output", help="Output directory (default: mirdp3_<timestamp>)")
    parser.add_argument("-t", "--temp", help="Temporary directory (default: output_dir/temp)")
    parser.add_argument("-clean", action="store_true", help="Remove temp and index after completion")
    args = parser.parse_args()

    # Step 1: check dependencies
    check_dependencies()

    # Step 2: determine base directory (where src/ and scripts/ are located)
    script_path = os.path.abspath(__file__)
    bin_dir = os.path.dirname(script_path)          # /home/bob/mirdp3/bin
    base_dir = os.path.dirname(bin_dir)             # /home/bob/mirdp3
    src_dir = os.path.join(base_dir, "src")
    scripts_dir = os.path.join(base_dir, "scripts")

    if not os.path.isdir(src_dir):
        sys.stderr.write(f"Error: src directory not found at {src_dir}\n")
        sys.exit(1)
    if not os.path.isdir(scripts_dir):
        sys.stderr.write(f"Error: scripts directory not found at {scripts_dir}\n")
        sys.exit(1)

    # Step 3: create output directory and temp directory
    if args.output:
        output_dir = args.output
    else:
        timestamp = datetime.now().strftime("%Y%m%d-%H%M")
        output_dir = f"mirdp3_{timestamp}"
    output_dir = os.path.abspath(ensure_dir(output_dir))

    # Set up log file
    log_file = os.path.join(output_dir, "mirdp3.log")
    with open(log_file, 'w') as lf:
        lf.write(f"Command: {full_command}\n")
        lf.write("=" * 80 + "\n")

    # Internal helper to run commands with logging
    def run_cmd(cmd, description):
        """Run a shell command, log output, print description to terminal."""
        print(f"[*] {description}")
        with open(log_file, 'a') as lf:
            lf.write(f"\n[{description}]\n")
            lf.write(f"Command: {cmd}\n")
            lf.write("Output:\n")
            lf.flush()
            # Use bash -c to support process substitution <( )
            full_cmd = f"bash -c {shlex.quote(cmd)}"
            ret = subprocess.call(full_cmd, shell=True, stdout=lf, stderr=subprocess.STDOUT)
            if ret != 0:
                lf.write(f"Error: command failed with exit code {ret}\n")
                sys.stderr.write(f"Error: command failed (see {log_file} for details)\n")
                sys.exit(1)
            lf.write("\n" + "-" * 40 + "\n")
        return ret

    if args.temp:
        temp_dir = args.temp
    else:
        temp_dir = os.path.join(output_dir, "temp")
    temp_dir = os.path.abspath(ensure_dir(temp_dir))

    # Index directory for *new* indexes (used in later steps)
    output_index_dir = os.path.join(output_dir, "index")
    ensure_dir(output_index_dir)

    # Step 4: determine existing index for PmiREN-core (if provided)
    pmiren_fa = os.path.abspath(args.pmiren)
    input_fa = os.path.abspath(args.input)

    if args.index:
        existing_index_prefix = os.path.abspath(args.index)
        if not check_index_exists(existing_index_prefix):
            sys.stderr.write(f"Error: Bowtie index does not exist at {existing_index_prefix}\n")
            sys.exit(1)
        use_existing_index = True
    else:
        existing_index_prefix = None
        use_existing_index = False

    # Define common file names
    aln_file = os.path.join(temp_dir, "isoform-out.aln")
    aln_fa = os.path.join(temp_dir, "isoform-out-aln.fa")
    nonaln_fa = os.path.join(temp_dir, "isoform-out-nonaln.fa")
    aln_length = os.path.join(temp_dir, "isoform-out-aln.length")
    pmiren_length = os.path.join(temp_dir, "isoform-in.length")
    temp_length = os.path.join(temp_dir, "isoform-temp.length")
    score_file = os.path.join(temp_dir, "isoform-out.score")
    best_score = os.path.join(temp_dir, "isoform-out.best.score")
    cluster_file = os.path.join(temp_dir, "isoform-out.best.score.cluster")

    if not use_existing_index:
        # Build index for PmiREN-core inside output_index_dir
        index_base = os.path.join(output_index_dir, "isoform-pmiren-core")
        run_cmd(f"bowtie-build -f {pmiren_fa} {index_base}",
                "Building Bowtie index for PmiREN-core dataset")
        bowtie_cmd = f"bowtie -a -v 2 {index_base} -f {input_fa} > {aln_file}"
    else:
        # Use existing index prefix
        bowtie_cmd = f"bowtie -a -v 2 {existing_index_prefix} -f {input_fa} > {aln_file}"
    run_cmd(bowtie_cmd, "Bowtie alignment of input against PmiREN-core")

    # Step 5: extract aligned and non-aligned sequences
    aligned_names = os.path.join(temp_dir, "aligned_names.txt")
    run_cmd(f"cat {aln_file} | awk '{{ print $1 }}' | sort -u > {aligned_names}",
            "Extracting aligned names list")
    run_cmd(f"seqkit grep -n -f {aligned_names} {input_fa} -o {aln_fa}",
            "Extracting aligned sequences")
    run_cmd(f"seqkit grep -n -v -f {aligned_names} {input_fa} -o {nonaln_fa}",
            "Extracting non-aligned sequences")

    # Step 6: process aligned part
    run_cmd(f"seqkit fx2tab --length --name --header-line {aln_fa} | grep -v '#' > {aln_length}",
            "Generating length table for aligned FASTA")
    run_cmd(f"seqkit fx2tab --length --name --header-line {pmiren_fa} | grep -v '#' > {pmiren_length}",
            "Generating length table for PmiREN-core")
    run_cmd(f"cat {aln_length} {pmiren_length} > {temp_length}",
            "Merging length tables")
    run_cmd(f"python {src_dir}/bowtie_mismatch_scoring.py -i {aln_file} -l {temp_length} -o {score_file}",
            "Computing mismatch scores for aligned reads")
    run_cmd(f"python {scripts_dir}/filter_best_scores_Q.py -i {score_file} -o {best_score}",
            "Filtering best scores")
    run_cmd(f"python {scripts_dir}/merge_by_scores_Q_v2.py -i {best_score} -o {cluster_file} -f {input_fa}",
            "Merging by scores (clustering)")

    # Step 7: process non-aligned part (first pass)
    nonaln_aln_file = os.path.join(temp_dir, "isoform-out-nonaln.aln")
    nonaln_length = os.path.join(temp_dir, "isoform-out-nonaln.length")
    nonaln_temp_length = os.path.join(temp_dir, "isoform-out-temp.length")
    nonaln_score = os.path.join(temp_dir, "isoform-out-nonaln.score")
    nonaln_best = os.path.join(temp_dir, "isoform-out-nonaln.best.score")

    nonaln_index_base = os.path.join(output_index_dir, "isoform-out")
    run_cmd(f"bowtie-build -f {nonaln_fa} {nonaln_index_base}",
            "Building Bowtie index for non-aligned sequences")
    run_cmd(f"bowtie -a -v 2 {nonaln_index_base} -f {nonaln_fa} > {nonaln_aln_file}",
            "Aligning non-aligned sequences against themselves")
    run_cmd(f"seqkit fx2tab --length --name --header-line {nonaln_fa} | grep -v '#' > {nonaln_length}",
            "Length table for non-aligned FASTA")
    run_cmd(f"cat {nonaln_length} {pmiren_length} > {nonaln_temp_length}",
            "Merging length tables (non-aligned + PmiREN)")
    run_cmd(f"python {src_dir}/bowtie_mismatch_scoring.py -i {nonaln_aln_file} -l {nonaln_temp_length} -o {nonaln_score}",
            "Scoring non-aligned alignments")
    run_cmd(f"python {scripts_dir}/filter_best_scores_Q.py -i {nonaln_score} -o {nonaln_best}",
            "Filtering best scores for non-aligned")

    # Step 8: split non-aligned into single and multiple best matches
    single_names_file = os.path.join(temp_dir, "isoform-out-nonaln.best.single")
    multi_names_file = os.path.join(temp_dir, "isoform-out-nonaln.best.multiple")
    single_fa = os.path.join(temp_dir, "isoform-out-nonaln-single.fa")
    multi_fa = os.path.join(temp_dir, "isoform-out-nonaln-multiple.fa")

    cmd_single = f"awk '{{ print $1 }}' {nonaln_best} | awk '{{sum[$1]+=1}}END{{for(i in sum) if(sum[i]==1) print i}}' > {single_names_file}"
    run_cmd(cmd_single, "Extracting sequences with single best hit (non-aligned)")
    run_cmd(f"seqkit grep -n -f {single_names_file} {nonaln_fa} -o {single_fa}",
            "Writing single-hit sequences")

    cmd_multi = f"awk '{{ print $1 }}' {nonaln_best} | awk '{{sum[$1]+=1}}END{{for(i in sum) if(sum[i]!=1) print i}}' > {multi_names_file}"
    run_cmd(cmd_multi, "Extracting sequences with multiple best hits")
    run_cmd(f"seqkit grep -n -f {multi_names_file} {nonaln_fa} -o {multi_fa}",
            "Writing multiple-hit sequences")

    # Step 9: process multiple-hit sequences against single-hit index
    single_index_base = os.path.join(output_index_dir, "isoform-out-single")
    run_cmd(f"bowtie-build -f {single_fa} {single_index_base}",
            "Building Bowtie index from single-hit sequences")
    multi_aln_file = os.path.join(temp_dir, "isoform-out-nonaln-multiple.aln")
    run_cmd(f"bowtie -a -v 2 {single_index_base} -f {multi_fa} > {multi_aln_file}",
            "Aligning multiple-hit sequences against single-hit index")

    multi_score = os.path.join(temp_dir, "isoform-out-nonaln-multiple.score")
    multi_best = os.path.join(temp_dir, "isoform-out-nonaln-multiple.best.score")
    single_length = os.path.join(temp_dir, "isoform-out-single.length")
    run_cmd(f"seqkit fx2tab --length --name --header-line {single_fa} | grep -v '#' > {single_length}",
            "Length table for single-hit FASTA")
    multi_temp_length = os.path.join(temp_dir, "isoform-out-multi-temp.length")
    run_cmd(f"cat {single_length} {nonaln_length} > {multi_temp_length}",
            "Merging length tables (single + multi)")
    run_cmd(f"python {src_dir}/bowtie_mismatch_scoring.py -i {multi_aln_file} -l {multi_temp_length} -o {multi_score}",
            "Scoring multiple-vs-single alignment")
    run_cmd(f"python {scripts_dir}/filter_best_scores_Q.py -i {multi_score} -o {multi_best}",
            "Filtering best scores for multiple-vs-single")

    # Extract sequences that are not aligned (become new singles)
    aligned_multi_names = os.path.join(temp_dir, "aligned_multi_names.txt")
    run_cmd(f"cat {multi_aln_file} | awk '{{ print $1 }}' | sort -u > {aligned_multi_names}",
            "Extracting aligned names from multiple-vs-single alignment")

    new_single_names = os.path.join(temp_dir, "isoform-out-nonaln-multiple.single")
    cmd_new_single = f"comm -23 <(cat {multi_fa} | grep '>' | sed 's/>//g' | sort) <(sort {aligned_multi_names}) > {new_single_names}"
    run_cmd(cmd_new_single, "Finding multi-hit sequences that did not align (become new singles)")

    total_single_names = os.path.join(temp_dir, "isoform-out-nonaln-total.single")
    run_cmd(f"cat {single_names_file} {new_single_names} > {total_single_names}",
            "Combining all single-hit sequence names")
    total_single_fa = os.path.join(temp_dir, "isoform-out-nonaln-total.single.fa")
    run_cmd(f"seqkit grep -n -f {total_single_names} {nonaln_fa} -o {total_single_fa}",
            "Extracting all single-hit sequences")

    multi_remaining_names = os.path.join(temp_dir, "isoform-out-nonaln-multiple.multiple")
    cmd_multi_remain = f"comm -23 <(cat {multi_fa} | grep '>' | sed 's/>//g' | sort) <(sort {new_single_names}) > {multi_remaining_names}"
    run_cmd(cmd_multi_remain, "Extracting sequences that remain multiple after alignment")
    multi_remaining_fa = os.path.join(temp_dir, "isoform-out-nonaln-multiple.multiple.fa")
    run_cmd(f"seqkit grep -n -f {multi_remaining_names} {multi_fa} -o {multi_remaining_fa}",
            "Writing remaining multiple-hit sequences")

    # Step 10: rename single sequences and align remaining multiple against them
    renamed_single_fa = os.path.join(temp_dir, "isoform-out-nonaln-total.single.rename.fa")
    run_cmd(f"python {src_dir}/rename_mir_fasta.py -i {total_single_fa} -s {args.start} -o {renamed_single_fa}",
            "Renaming single-hit sequences with start number")

    renamed_index_base = os.path.join(output_index_dir, "isoform-out-single-rename")
    run_cmd(f"bowtie-build -f {renamed_single_fa} {renamed_index_base}",
            "Building Bowtie index from renamed single sequences")
    multi_vs_rename_aln = os.path.join(temp_dir, "isoform-out-nonaln-multiple-single.aln")
    run_cmd(f"bowtie -a -v 2 {renamed_index_base} -f {multi_remaining_fa} > {multi_vs_rename_aln}",
            "Aligning remaining multiple sequences against renamed single index")

    rename_single_length = os.path.join(temp_dir, "isoform-single-rename.length")
    run_cmd(f"seqkit fx2tab --length --name --header-line {renamed_single_fa} | grep -v '#' > {rename_single_length}",
            "Length table for renamed single FASTA")
    multi_remain_length = os.path.join(temp_dir, "isoform-out-multi-remain.length")
    run_cmd(f"seqkit fx2tab --length --name --header-line {multi_remaining_fa} | grep -v '#' > {multi_remain_length}",
            "Length table for remaining multiple sequences")
    rename_total_length = os.path.join(temp_dir, "isoform-rename-total.length")
    run_cmd(f"cat {rename_single_length} {multi_remain_length} > {rename_total_length}",
            "Merging length tables for rename alignment")
    multi_rename_score = os.path.join(temp_dir, "isoform-out-nonaln-multiple-single.score")
    multi_rename_best = os.path.join(temp_dir, "isoform-out-nonaln-multiple-single.best.score")
    multi_rename_cluster = os.path.join(temp_dir, "isoform-out-nonaln-multiple-single.best.score.cluster")
    run_cmd(f"python {src_dir}/bowtie_mismatch_scoring.py -i {multi_vs_rename_aln} -l {rename_total_length} -o {multi_rename_score}",
            "Scoring alignment of remaining multiple vs renamed single")
    run_cmd(f"python {scripts_dir}/filter_best_scores_Q.py -i {multi_rename_score} -o {multi_rename_best}",
            "Filtering best scores")
    run_cmd(f"python {scripts_dir}/merge_by_scores_Q_v2.py -i {multi_rename_best} -o {multi_rename_cluster} -f {input_fa}",
            "Merging by scores for rename alignment")

    # Step 11: rename final sequences using mapping files
    renamed_aln_fa = os.path.join(temp_dir, "isoform-out-aln-rename.fa")
    run_cmd(f"python {src_dir}/rename_fasta_by_mapping.py -m {cluster_file} -i {aln_fa} -o {renamed_aln_fa}",
            "Renaming aligned sequences based on cluster mapping")
    renamed_multi_fa = os.path.join(temp_dir, "isoform-out-nonaln-multiple.multiple.renamed.fa")
    run_cmd(f"python {src_dir}/rename_fasta_by_mapping.py -m {multi_rename_cluster} -i {multi_remaining_fa} -o {renamed_multi_fa}",
            "Renaming multiple sequences based on cluster mapping")

    # Step 12: combine final outputs
    final_renamed = os.path.join(output_dir, "isoform-out-total.renamed.fa")
    final_new = os.path.join(output_dir, "isoform-out-total-new.renamed.fa")
    run_cmd(f"cat {renamed_aln_fa} {renamed_multi_fa} > {final_renamed}",
            "Combining renamed aligned and renamed multiple sequences")
    run_cmd(f"cat {renamed_single_fa} > {final_new}",
            "Copying renamed single sequences as final new sequences")

    print(f"\n[SUCCESS] Pipeline completed. Results saved in {output_dir}")
    print(f"  - Combined renamed sequences: {final_renamed}")
    print(f"  - New renamed single sequences: {final_new}")
    print(f"  - Log file: {log_file}")

    if args.clean:
        print("[*] Cleaning temporary and index directories...")
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
            print(f"  Removed {temp_dir}")
        if not args.index and os.path.exists(output_index_dir):
            shutil.rmtree(output_index_dir)
            print(f"  Removed {output_index_dir}")
        print("[*] Cleanup complete.")

if __name__ == "__main__":
    main()
