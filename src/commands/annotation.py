#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
MirDeep-P3 annotation step orchestrator.
Handles input (identification output folders), dependency checks,
output directory structure, replicate grouping, and the complete
annotation pipeline (4.4 - 4.13), with proper consistency handling.
"""

import argparse
import os
import sys
import shutil
import subprocess
import multiprocessing
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional, Tuple

from utils.config import load_config
from utils.dependencies import check_external_tools
from check_python_deps import check_python_deps

DEFAULT_THREADS = 1
DEFAULT_PROGRESS = 1

# ----------------------------------------------------------------------
def add_arguments(parser: argparse.ArgumentParser):
    """Define all command-line arguments for the annotation step."""
    parser.add_argument("-i", "--input", help="Input folder(s) from identification step, comma separated")
    parser.add_argument("-f", "--file", help="File containing list of input folders, one per line")   # 新增
    parser.add_argument("-o", "--output", help="Output root directory")
    parser.add_argument("--prefix", help="Output prefix(es) comma separated")
    parser.add_argument("--prefix_miRNA", help="miRNA prefix (e.g., 'N710') [required]")
    parser.add_argument("--species", help="Species name (e.g., 'Arabidopsis thaliana') [required]")

    parser.add_argument("-g", "--genome", help="Reference genome FASTA file")
    parser.add_argument("-d", "--index", help="Bowtie genome index prefix (skip building if provided)")

    parser.add_argument("-r", "--replicate", help="Replicate grouping: comma-separated counts per group")

    parser.add_argument("-t", "--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("-p", "--progress", type=int, default=DEFAULT_PROGRESS,
                        help="Number of parallel processes")

    parser.add_argument("-s", "--same", action="store_true",
                        help="Use same prefix for input and output (output folders will be 'case1','case2',...)")
    parser.add_argument("--consistency", help="Consistency reference basic-info file (e.g., YL1-basic-info)")
    parser.add_argument("--common", action="store_true", help="Enable common miRNA consistency analysis")
    parser.add_argument("--clean", action="store_true", help="Remove temporary directories (temp, index) after completion")

    parser.add_argument("--bowtie", help="path to bowtie")
    parser.add_argument("--bowtie-build", help="path to bowtie-build")
    parser.add_argument("--RNAfold", help="path to RNAfold")
    parser.add_argument("--samtools", help="path to samtools")
    parser.add_argument("--bedtools", help="path to bedtools")
    parser.add_argument("--seqkit", help="path to seqkit")

    parser.add_argument("-h", "--help", action="store_true", help="Show this help message and exit.")


def ensure_defaults(args):
    """Ensure all annotation-specific attributes exist on args with safe defaults."""
    defaults = {
        'prefix_miRNA': None,
        'species': None,
        'genome': None,
        'index': None,
        'replicate': '1',
        'threads': 1,
        'progress': 1,
        'prefix': None,
        'same': False,
        'consistency': None,
        'common': False,
        'clean': False,
        'bowtie': None,
        'bowtie_build': None,
        'RNAfold': None,
        'samtools': None,
        'bedtools': None,
        'seqkit': None,
        'input': None,
        'output': None,
        'file': None,                   # 新增
    }
    for attr, val in defaults.items():
        if not hasattr(args, attr):
            setattr(args, attr, val)


def validate_input_folders(folders: List[Path]):
    """Verify each input folder exists and contains required identification output files."""
    for folder in folders:
        if not folder.is_dir():
            sys.exit(f"Input folder not found: {folder}")
        pred_files = list(folder.glob("*_filter_P_prediction"))
        if not pred_files:
            sys.exit(f"No *_filter_P_prediction found in {folder}. Ensure identification step completed.")
    print("Input folder validation passed.")


def detect_file_prefix(input_folder: Path) -> str:
    """Extract the file prefix from the identification output folder."""
    pred_files = list(input_folder.glob("*_filter_P_prediction"))
    if pred_files:
        return pred_files[0].stem.replace("_filter_P_prediction", "")
    # Fallback: folder name
    return input_folder.name


def resolve_annotation_prefixes(
    replicate_groups: List[List[Path]],
    user_prefix: Optional[str],
    same: bool,
    input_folders: List[Path]
) -> Tuple[List[str], List[str], List[Path]]:
    """
    Compute output folder prefixes, file prefixes, and output directories per input.
    Returns:
        folder_prefixes: list of folder names (one per group)
        file_prefixes:   list of file prefixes (one per input)
        output_dirs:     list of output paths (one per input)
    """
    num_groups = len(replicate_groups)
    num_inputs = len(input_folders)

    # ---- 1. Determine file prefixes ----
    if same:
        # same mode: user must provide one prefix per input
        if not user_prefix:
            sys.exit("Error: --same requires --prefix with one prefix per input file.")
        file_prefixes = [p.strip() for p in user_prefix.split(",")]
        if len(file_prefixes) != num_inputs:
            sys.exit(f"Error: --same requires --prefix to have exactly {num_inputs} entries (got {len(file_prefixes)}).")
    else:
        if user_prefix is None:
            # Auto-detect file prefixes for each input folder
            file_prefixes = [detect_file_prefix(f) for f in input_folders]
        else:
            prefix_candidates = [p.strip() for p in user_prefix.split(",")]
            if len(prefix_candidates) == 1:
                # 1 prefix -> it becomes the output folder name; file prefixes auto-detected
                file_prefixes = [detect_file_prefix(f) for f in input_folders]
            elif len(prefix_candidates) == num_inputs:
                # One prefix per input -> these are the file prefixes
                file_prefixes = prefix_candidates
            elif len(prefix_candidates) == num_groups:
                # One prefix per group -> each group gets its own folder, file prefixes auto-detected
                file_prefixes = [detect_file_prefix(f) for f in input_folders]
            else:
                sys.exit(f"Error: --prefix count ({len(prefix_candidates)}) must be 1, {num_groups} (groups), or {num_inputs} (inputs).")

    # ---- 2. Determine folder prefixes (output directory names) ----
    if same:
        folder_prefixes = [f"case{i+1}" for i in range(num_groups)]
    else:
        if user_prefix is None:
            # No user prefix -> auto-generate case1, case2...
            folder_prefixes = [f"case{i+1}" for i in range(num_groups)]
        else:
            prefix_candidates = [p.strip() for p in user_prefix.split(",")]
            if len(prefix_candidates) == 1:
                # Single prefix -> all groups use this same folder
                base = prefix_candidates[0]
                folder_prefixes = [base] * num_groups
            elif len(prefix_candidates) == num_groups:
                # One folder name per group
                folder_prefixes = prefix_candidates
            elif len(prefix_candidates) == num_inputs:
                # File prefixes specified individually -> folders are case1, case2...
                folder_prefixes = [f"case{i+1}" for i in range(num_groups)]
            else:
                # Should not reach here
                sys.exit("Unexpected prefix count.")

    # ---- 3. Map each input to its output directory ----
    output_dirs = []
    for g_idx, grp in enumerate(replicate_groups):
        for _ in grp:
            output_dirs.append(Path(folder_prefixes[g_idx]))

    # Final consistency check
    assert len(file_prefixes) == num_inputs
    assert len(output_dirs) == num_inputs
    return folder_prefixes, file_prefixes, output_dirs


# ----------------------------------------------------------------------
# Part 1: steps 4.4 – 4.10, returns path to second-basic-info file
def run_annotation_part1(
    args,
    group_input_folders: List[Path],
    group_file_prefixes: List[str],
    output_group_dir: Path,
    genome: Path,
    log_file: Path,
    project_root: Path
) -> Path:
    """
    Execute annotation steps 4.4 to 4.10 for one replicate group.
    Returns the path to {prefix}-second-basic-info file (in temp/).
    """
    py_exe = sys.executable
    src_dir = project_root / "src"
    scripts_dir = project_root / "scripts"
    data_dir = project_root / "data"
    temp_dir = output_group_dir / "temp"
    index_dir = output_group_dir / "index"
    temp_dir.mkdir(parents=True, exist_ok=True)
    index_dir.mkdir(parents=True, exist_ok=True)

    num_samples = len(group_input_folders)
    single = (num_samples == 1)
    # use the first file prefix as the group prefix for combined outputs
    group_prefix = output_group_dir.name  # e.g., "case1", "rootAA"

    # 4.4 Calculate RPM and generate all-mod predictions
    read_total_file = output_group_dir / "read_total_file"
    with open(read_total_file, 'w') as rtf:
        for folder, fprefix in zip(group_input_folders, group_file_prefixes):
            total_read_path = folder / f"{fprefix}.total_reads"
            if not total_read_path.exists():
                sys.exit(f"Missing total_reads file: {total_read_path}")
            with open(total_read_path) as f:
                count = f.read().strip()
            rtf.write(f"{fprefix}\t{count}\n")

    all_mod_files = []
    for folder, fprefix in zip(group_input_folders, group_file_prefixes):
        pred_input = folder / f"{fprefix}_filter_P_prediction"
        out_all_pred = output_group_dir / f"{fprefix}-all-mod_fp_prediction"
        count_to_rpm = temp_dir / f"{fprefix}_count_to_rpm"
        subprocess.run([py_exe, str(src_dir / "mirdp3_all_pred.py"),
                        "-i", str(pred_input), "-c", str(read_total_file),
                        "-o", str(out_all_pred), "-m", str(count_to_rpm)],
                       stdout=subprocess.DEVNULL, stderr=open(log_file, 'a'), check=True)
        all_mod_files.append(out_all_pred)

    if single:
        merged_pred = all_mod_files[0]
    else:
        merged_pred = output_group_dir / f"{group_prefix}-all-mod_fp_prediction"
        with open(merged_pred, 'w') as mf:
            for fpath in all_mod_files:
                with open(fpath) as inf:
                    mf.write(inf.read())

    # 4.5 Convert to BED, merge, remove redundancy, extract mature
    bed_tmp = temp_dir / f"{group_prefix}-all-mod_fp_prediction-tmp.bed"
    subprocess.run([py_exe, str(src_dir / "all_pred_to_bed.py"),
                    "-i", str(merged_pred), "-o", str(bed_tmp)],
                   stdout=subprocess.DEVNULL, stderr=open(log_file, 'a'), check=True)

    sorted_bed = temp_dir / f"{group_prefix}-all-mod_fp_prediction-tmp.sorted.bed"
    subprocess.run([py_exe, str(src_dir / "all_pred_to_bed_sort.py"),
                    "-i", str(bed_tmp), "-o", str(sorted_bed)],
                   stdout=subprocess.DEVNULL, stderr=open(log_file, 'a'), check=True)

    merged_bed = temp_dir / f"{group_prefix}-merge_output.bed"
    bedtools_exe = getattr(args, 'bedtools', None) or "bedtools"
    subprocess.run(f"{bedtools_exe} merge -d 0 -c 4,5,6 -o collapse,distinct,distinct -i {sorted_bed} > {merged_bed}",
                   shell=True, check=True, stdout=subprocess.DEVNULL, stderr=open(log_file, 'a'))

    nr_pred = temp_dir / f"{group_prefix}-all-mod_fp_prediction-nr"
    subprocess.run([py_exe, str(src_dir / "remove_nr_reads.py"),
                    "-m", str(merged_bed), "-i", str(merged_pred), "-o", str(nr_pred)],
                   stdout=subprocess.DEVNULL, stderr=open(log_file, 'a'), check=True)

    mature_bed = temp_dir / f"{group_prefix}-all-mod_fp_prediction-nr-mature.bed"
    mature_fasta = temp_dir / f"{group_prefix}-all-mod_fp_prediction-nr-mature.fasta"
    subprocess.run([py_exe, str(src_dir / "extract_mature.py"),
                    "-i", str(nr_pred), "-genome", str(genome),
                    "-o", str(mature_bed), "-fo", str(mature_fasta)],
                   stdout=subprocess.DEVNULL, stderr=open(log_file, 'a'), check=True)

    # 4.6 Align to PmiREN
    pmiren_index = data_dir / "index" / "PmiREN-core-AA"
    aln_file = temp_dir / f"{group_prefix}-all-mod_fp_prediction-nr-mature-PmiREN.aln"
    bowtie_exe = getattr(args, 'bowtie', None) or "bowtie"
    subprocess.run(f"{bowtie_exe} -a -v 2 {pmiren_index} -f {mature_fasta} > {aln_file} 2>> {log_file}",
                   shell=True, check=True)

    # 4.7 Process alignment, scoring, clustering, annotate
    seqkit_exe = getattr(args, 'seqkit', None) or "seqkit"
    len_file = temp_dir / f"{group_prefix}-all-mod_fp_prediction-nr-mature.length"
    subprocess.run(f"{seqkit_exe} fx2tab --length --name --header-line {mature_fasta} | grep -v '#' > {len_file} 2>> {log_file}",
                   shell=True, check=True)

    pmiren_len_file = data_dir / "PmiREN-core-in-20260429-AA.length"
    total_len = temp_dir / f"{group_prefix}-all-mod_fp_prediction-nr-mature-total.length"
    with open(total_len, 'w') as outf:
        for f in [pmiren_len_file, len_file]:
            with open(f) as inf:
                outf.write(inf.read())

    score_file = temp_dir / f"{group_prefix}-all-mod_fp_prediction-nr-mature-PmiREN.score"
    subprocess.run([py_exe, str(src_dir / "bowtie_mismatch_scoring.py"),
                    "-i", str(aln_file), "-l", str(total_len), "-o", str(score_file)],
                   stdout=subprocess.DEVNULL, stderr=open(log_file, 'a'), check=True)

    best_score = temp_dir / f"{group_prefix}-all-mod_fp_prediction-nr-mature-PmiREN.best.score"
    subprocess.run([py_exe, str(scripts_dir / "filter_best_scores_Q.py"),
                    "-i", str(score_file), "-o", str(best_score)],
                   stdout=subprocess.DEVNULL, stderr=open(log_file, 'a'), check=True)

    cluster_file = temp_dir / f"{group_prefix}-all-mod_fp_prediction-nr-mature-PmiREN.best.score.cluster"
    subprocess.run([py_exe, str(scripts_dir / "merge_by_scores_Q_v2.py"),
                    "-i", str(best_score), "-o", str(cluster_file), "-f", str(mature_fasta)],
                   stdout=subprocess.DEVNULL, stderr=open(log_file, 'a'), check=True)

    annotated_nr = temp_dir / f"{group_prefix}-all-mod_fp_prediction-nr-annotated"
    subprocess.run([py_exe, str(src_dir / "predictions_annotate.py"),
                    "-i", str(nr_pred), "-c", str(cluster_file), "-o", str(annotated_nr)],
                   stdout=subprocess.DEVNULL, stderr=open(log_file, 'a'), check=True)

    # 4.8 Filter non-conserved
    filtered_nr = temp_dir / f"{group_prefix}-filtered-all-nr-anno"
    subprocess.run([py_exe, str(src_dir / "predictions_annotate_filter.py"),
                    "-i", str(annotated_nr), "-o", str(filtered_nr)],
                   stdout=subprocess.DEVNULL, stderr=open(log_file, 'a'), check=True)

    # 4.9 Extract structures
    struct_file = temp_dir / f"{group_prefix}-stemloop.struc"
    struct_20nt = temp_dir / f"{group_prefix}-stemloop_20nt.struc"
    fasta_20nt = temp_dir / f"{group_prefix}-stemloop_20nt.fasta"
    subprocess.run([py_exe, str(src_dir / "extract_struc.py"),
                    "-i", str(filtered_nr), "-genome", str(genome),
                    "--threads", str(args.threads),
                    "-struc", str(struct_file), "-struc_20nt", str(struct_20nt),
                    "-fasta_20nt", str(fasta_20nt)],
                   stdout=subprocess.DEVNULL, stderr=open(log_file, 'a'), check=True)

    # 4.10 Primary basic info + rearrange → second-basic-info
    chr_len_file = group_input_folders[0] / "chr_length"
    if not chr_len_file.exists():
        sys.exit(f"chr_length file missing in {group_input_folders[0]}")
    primary_info = temp_dir / f"{group_prefix}-primary-basic-info"
    subprocess.run([py_exe, str(src_dir / "primary-basic-info.py"),
                    "-f", str(filtered_nr), "-struc_20nt", str(struct_20nt),
                    "-struc", str(struct_file), "-chr_length", str(chr_len_file),
                    "-species", args.species,
                    "-prefix_miRNA", args.prefix_miRNA,
                    "-o", str(primary_info)],
                   stdout=subprocess.DEVNULL, stderr=open(log_file, 'a'), check=True)

    second_basic_info = temp_dir / f"{group_prefix}-second-basic-info"
    subprocess.run([py_exe, str(src_dir / "rearrange_result.py"),
                    "-i", str(primary_info), "-o", str(second_basic_info)],
                   stdout=subprocess.DEVNULL, stderr=open(log_file, 'a'), check=True)

    return second_basic_info  # temp/{prefix}-second-basic-info

# ----------------------------------------------------------------------
# MIRN rename step (between 4.10 and 4.11)
def run_mirn_rename(
    args,
    second_basic_info: Path,
    output_group_dir: Path,
    prefix: str,
    log_file: Path,
    project_root: Path
) -> Path:
    """
    Execute MIRN rename pipeline.
    Uses MIRN_cluster.py, MIRN_rename.py, add_AA.py, etc.
    Returns path to final basic-info file (output/{prefix}-basic-info).
    """
    py_exe = sys.executable
    src_dir = project_root / "src"
    scripts_dir = project_root / "scripts"
    bin_dir = project_root / "bin"
    temp_dir = output_group_dir / "temp"
    index_dir = output_group_dir / "index"
    temp_dir.mkdir(parents=True, exist_ok=True)
    index_dir.mkdir(parents=True, exist_ok=True)

    # Prepare MIRN sequences
    mirn_fa = temp_dir / f"{prefix}-second-basic-info-MIRN.fa"
    subprocess.run(f"cat {second_basic_info} | grep \"MIRN\" | awk '{{ print \">\"$1\"\\n\"$18 }}' > {mirn_fa}",
                   shell=True, check=True, stderr=open(log_file, 'a'))

    mirn_len = temp_dir / f"{prefix}-second-basic-info-MIRN.length"
    seqkit_exe = getattr(args, 'seqkit', None) or "seqkit"
    subprocess.run(f"{seqkit_exe} fx2tab --length --name --header-line {mirn_fa} | grep -v '#' > {mirn_len} 2>> {log_file}",
                   shell=True, check=True)

    mirn_a_fa = temp_dir / f"{prefix}-second-basic-info-MIRN-A.fa"
    subprocess.run([py_exe, str(bin_dir / "add_AA.py"),
                    "-i", str(mirn_fa), "-o", str(mirn_a_fa)],
                   stdout=subprocess.DEVNULL, stderr=open(log_file, 'a'), check=True)

    # Build bowtie index
    bowtie_build_exe = getattr(args, 'bowtie_build', None) or "bowtie-build"
    mirn_index = index_dir / "basic-info-MIRN-A"
    subprocess.run(f"{bowtie_build_exe} -f {mirn_a_fa} {mirn_index} >> {log_file} 2>&1",
                   shell=True, check=True)

    # Align original MIRN to AA extended
    bowtie_exe = getattr(args, 'bowtie', None) or "bowtie"
    mirn_aln = temp_dir / f"{prefix}-second-basic-info-MIRN.aln"
    subprocess.run(f"{bowtie_exe} -a -v 2 {mirn_index} -f {mirn_fa} > {mirn_aln} 2>> {log_file}",
                   shell=True, check=True)

    # Find single-copy MIRN
    mirn_single = temp_dir / f"{prefix}-second-basic-info-MIRN.single"
    subprocess.run(f"cat {mirn_aln} | awk '{{sum[$1]+=1}}END{{for(i in sum){{ print i\"\\t\"sum[i] }}}}' | awk '$2==1{{ print $1 }}' > {mirn_single}",
                   shell=True, check=True, stderr=open(log_file, 'a'))

    # MIRN that align to different sequences (not self)
    mirn_single_part = temp_dir / f"{prefix}-second-basic-info-MIRN.single.part"
    subprocess.run(f"cat {mirn_aln} | awk '{{sum[$1]+=1}}END{{for(i in sum){{ print i\"\\t\"sum[i] }}}}' | awk '$2!=1{{ print $1 }}' | awk 'NR==FNR {{a[$1];next}} {{if($1 in a){{ print $0 }}}}' - {mirn_aln} | awk '$1!=$3{{ print $0 }}' | awk 'NR==FNR {{a[$1];next}} {{if($3 in a){{ print $0 }}}}' {mirn_single} - | awk '{{ print $1\"\\t\"$3 }}' > {mirn_single_part}",
                   shell=True, check=True, stderr=open(log_file, 'a'))

    # MIRN alignments excluding self-alignments (for scoring)
    mirn_aln_part = temp_dir / f"{prefix}-second-basic-info-MIRN.aln.part"
    subprocess.run(f"cat {mirn_aln} | awk '{{sum[$1]+=1}}END{{for(i in sum){{ print i\"\\t\"sum[i] }}}}' | awk '$2!=1{{ print $1 }}' | awk 'NR==FNR {{a[$1];next}} {{if($1 in a){{ print $0 }}}}' - {mirn_aln} | awk '$1!=$3{{ print $0 }}' | awk 'NR==FNR {{a[$1];next}} {{if(!($3 in a)){{ print $0 }}}}' {mirn_single} - > {mirn_aln_part}",
                   shell=True, check=True, stderr=open(log_file, 'a'))

    # Score, filter best, cluster
    score_file = temp_dir / f"{prefix}-second-basic-info-MIRN.aln.part.score"
    subprocess.run([py_exe, str(src_dir / "bowtie_mismatch_scoring.py"),
                    "-i", str(mirn_aln_part), "-l", str(mirn_len), "-o", str(score_file)],
                   stdout=subprocess.DEVNULL, stderr=open(log_file, 'a'), check=True)

    best_score = temp_dir / f"{prefix}-second-basic-info-MIRN.aln.part.best.score"
    subprocess.run([py_exe, str(scripts_dir / "filter_best_scores_Q.py"),
                    "-i", str(score_file), "-o", str(best_score)],
                   stdout=subprocess.DEVNULL, stderr=open(log_file, 'a'), check=True)

    cluster = temp_dir / f"{prefix}-second-basic-info-MIRN.aln.part.best.score.cl"
    subprocess.run([py_exe, str(src_dir / "MIRN_cluster.py"),
                    "-i", str(best_score), "-o", str(cluster)],
                   stdout=subprocess.DEVNULL, stderr=open(log_file, 'a'), check=True)

    # MIRN rename
    mirn_map = temp_dir / f"{prefix}-second-basic-info-MIRN.map"
    subprocess.run([py_exe, str(src_dir / "MIRN_rename.py"),
                    "-i", str(cluster), "-s", str(mirn_single),
                    "-p", str(mirn_single_part), "-o", str(mirn_map)],
                   stdout=subprocess.DEVNULL, stderr=open(log_file, 'a'), check=True)

    # Apply rename to second-basic-info → output/{prefix}-basic-info
    final_basic_info = output_group_dir / f"{prefix}-basic-info"
    subprocess.run([py_exe, str(src_dir / "rename-basic-info-MIRN.py"),
                    "-i", str(second_basic_info), "-m", str(mirn_map),
                    "-o", str(final_basic_info)],
                   stdout=subprocess.DEVNULL, stderr=open(log_file, 'a'), check=True)

    return final_basic_info

# Part 2: steps 4.12 – 4.13, using the given basic_info file
def run_annotation_part2(
    args,
    group_input_folders: List[Path],
    group_file_prefixes: List[str],
    output_group_dir: Path,
    genome: Path,
    basic_info_file: Path,
    project_root: Path
):
    """Execute 4.12 (cluster) and 4.13 (expression) using final basic-info."""
    py_exe = sys.executable
    src_dir = project_root / "src"
    temp_dir = output_group_dir / "temp"
    index_dir = output_group_dir / "index"
    temp_dir.mkdir(parents=True, exist_ok=True)

    # 4.12 miRNA cluster
    mature_fa_cluster = temp_dir / f"{output_group_dir.name}-basic-info-mature.fa"
    with open(basic_info_file) as fin, open(mature_fa_cluster, 'w') as fout:
        for line in fin:
            parts = line.strip().split('\t')
            if len(parts) >= 18:
                fout.write(f">{parts[0]}\n{parts[14]}\n{parts[17]}\n")

    out_cluster = output_group_dir / f"{output_group_dir.name}-basic-info-cluster"
    subprocess.run([py_exe, str(src_dir / "scan_cluster.py"),
                    "-i", str(basic_info_file), "-species", args.species,
                    "-prefix", args.prefix_miRNA, "-o", str(out_cluster)],
                   stdout=subprocess.DEVNULL, stderr=open(output_group_dir / f"{output_group_dir.name}_annotation.log", 'a'), check=True)

    # 4.13 Expression matrix
    mature_1nt_fa = temp_dir / f"{output_group_dir.name}-basic-info-mature-1nt.fa"
    subprocess.run([py_exe, str(src_dir / "extract_mature_flank.py"),
                    "-i", str(basic_info_file), "-o", str(mature_1nt_fa)],
                   stdout=subprocess.DEVNULL, stderr=open(output_group_dir / f"{output_group_dir.name}_annotation.log", 'a'), check=True)

    mature_index = index_dir / f"{output_group_dir.name}-mature-1nt"
    bowtie_build_exe = getattr(args, 'bowtie_build', None) or "bowtie-build"
    subprocess.run([bowtie_build_exe, "-f", str(mature_1nt_fa), str(mature_index)],
                   stdout=subprocess.DEVNULL, stderr=open(output_group_dir / f"{output_group_dir.name}_annotation.log", 'a'), check=True)

    bowtie_exe = getattr(args, 'bowtie', None) or "bowtie"
    aln_files = []
    for folder, fprefix in zip(group_input_folders, group_file_prefixes):
        candidates = list(folder.rglob("*.processed.fa"))
        if not candidates:
            sys.exit(f"No processed.fa found in {folder} or its subdirectories.")
        proc_fa = candidates[0]
        out_aln = temp_dir / f"{fprefix}.mature.aln"
        subprocess.run(f"{bowtie_exe} -a -v 0 {mature_index} -f {proc_fa} > {out_aln} 2>> {output_group_dir / f'{output_group_dir.name}_annotation.log'}",
                       shell=True, check=True)
        aln_files.append(out_aln)

    read_total_file = output_group_dir / "read_total_file"
    aln_list = ",".join(str(a) for a in aln_files)
    exp_prefix = output_group_dir.name
    subprocess.run([py_exe, str(src_dir / "expression_matrix.py"),
                    "-aln", aln_list, "-c", str(read_total_file),
                    "-species", args.species,
                    "-o", str(output_group_dir / f"{exp_prefix}-mature.exp"),
                    "-oc", str(output_group_dir / f"{exp_prefix}-mature.count")],
                   stdout=subprocess.DEVNULL, stderr=open(output_group_dir / f"{output_group_dir.name}_annotation.log", 'a'), check=True)


# ----------------------------------------------------------------------
def run(args):
    """Main entry point for annotation subcommand."""
    ensure_defaults(args)
    cfg = getattr(args, 'config_data', {})
    project_root = args.project_root

    # Merge config file values (same as identification)
    config_mapping = {
        'i/input': 'input',
        'o/output': 'output',
        'g/genome': 'genome',
        't/threads': 'threads',
        'p/progress': 'progress',
        'r/replicate': 'replicate',
        'd/index': 'index',
        'prefix': 'prefix',
        'prefix_miRNA': 'prefix_miRNA',
        'species': 'species',
        'bowtie': 'bowtie',
        'bowtie-build': 'bowtie_build',
        'RNAfold': 'RNAfold',
        'samtools': 'samtools',
        'bedtools': 'bedtools',
        'seqkit': 'seqkit',
        'consistency': 'consistency',
        'common': 'common',
        'clean': 'clean',
        'same': 'same',
        'f/file': 'file',
    }
    for cfg_key, attr in config_mapping.items():
        if cfg_key in cfg and getattr(args, attr, None) is None:
            setattr(args, attr, cfg[cfg_key])

    # Required argument check
    if not args.prefix_miRNA:
        sys.exit("Error: --prefix_miRNA is required.")
    if not args.species:
        sys.exit("Error: --species is required.")

    # Dependency checks
    required_tools = {
        'bowtie': getattr(args, 'bowtie', None) or cfg.get('bowtie'),
        'bowtie-build': getattr(args, 'bowtie_build', None) or cfg.get('bowtie-build'),
        'RNAfold': getattr(args, 'RNAfold', None) or cfg.get('RNAfold'),
        'samtools': getattr(args, 'samtools', None) or cfg.get('samtools'),
        'bedtools': getattr(args, 'bedtools', None) or cfg.get('bedtools'),
        'seqkit': getattr(args, 'seqkit', None) or cfg.get('seqkit'),
    }
    all_ok, missing = check_external_tools(required_tools)
    if not all_ok:
        sys.exit("Error: missing required external dependencies.")
    check_python_deps()

    # Input folders (support -i, -f, config)
    input_folders_raw = []
    if args.file:
        file_path = Path(args.file)
        if not file_path.is_file():
            sys.exit(f"Input list file not found: {args.file}")
        with open(file_path) as f:
            input_folders_raw = [line.strip() for line in f if line.strip()]
    elif args.input:
        input_folders_raw = [p.strip() for p in args.input.split(",")]
    else:
        raw_cfg = cfg.get("i/input") or cfg.get("input")
        if not raw_cfg:
            sys.exit("Error: input folders must be specified via -i, -f, or config.")
        if Path(raw_cfg).is_file():
            with open(raw_cfg) as f:
                input_folders_raw = [line.strip() for line in f if line.strip()]
        else:
            input_folders_raw = [p.strip() for p in raw_cfg.split(",")]
    input_folders = [Path(p) for p in input_folders_raw if p]
    validate_input_folders(input_folders)

    # Genome
    genome = Path(args.genome or cfg.get("g/genome") or cfg.get("genome"))
    if not genome.exists():
        sys.exit(f"Genome file not found: {genome}")

    # Replicate grouping
    rep_raw = args.replicate or cfg.get("r/replicate") or cfg.get("replicate") or "1"
    rep_counts = [int(x.strip()) for x in rep_raw.split(",")]
    replicate_groups = []
    idx = 0
    for cnt in rep_counts:
        replicate_groups.append(input_folders[idx:idx+cnt])
        idx += cnt
    if idx != len(input_folders):
        sys.exit("Error: replicate counts do not sum to number of input folders.")
    args.replicate_groups = replicate_groups
    num_groups = len(replicate_groups)

    # Output directory
    output_root = Path(args.output or cfg.get("o/output") or cfg.get("output") or f"mirdeep-annotation-{datetime.now().strftime('%m%d%y-%H%M')}")
    output_root.mkdir(parents=True, exist_ok=True)

    # Prefix resolution (auto-detect if not given)
    user_prefix = args.prefix or cfg.get("prefix")
    folder_prefixes, file_prefixes, output_dirs = resolve_annotation_prefixes(
        replicate_groups, user_prefix, args.same, input_folders
    )
    # Override file_prefixes if auto-detection returned None (should not happen)
    if file_prefixes is None:
        file_prefixes = [detect_file_prefix(f) for f in input_folders]

    # Pipe file with timestamp
    timestamp = datetime.now().strftime("%m%d%Y-%H%M")
    pipe_file = output_root / f"mirdp3-annotation-{timestamp}.pipe"
    with open(pipe_file, 'w') as pf:
        for inp, out_rel, fpref in zip(input_folders, output_dirs, file_prefixes):
            out_path = output_root / out_rel
            pf.write(f"{inp}\t{out_path}\t{fpref}\n")

    print(f"\nProcessing {len(input_folders)} input folder(s) with {args.progress} parallel process(es).")

    # ---- Execute part 1 (per group) ----
    second_basic_info_paths = []
    for g_idx, grp in enumerate(replicate_groups):
        out_dir = output_root / folder_prefixes[g_idx]
        group_fprefixes = [file_prefixes[input_folders.index(inp)] for inp in grp]
        log_file = out_dir / f"{folder_prefixes[g_idx]}_annotation.log"
        log_file.parent.mkdir(parents=True, exist_ok=True)

        print(f"  [start] group {g_idx+1} ({folder_prefixes[g_idx]})")
        seconds = run_annotation_part1(args, grp, group_fprefixes, out_dir, genome, log_file, project_root)
        second_basic_info_paths.append(seconds)
        print(f"  [done]  group {g_idx+1} ({folder_prefixes[g_idx]})")

    # ---- MIRN rename for each group ----
    basic_info_paths = []
    for g_idx, (out_dir, seconds) in enumerate(zip([output_root / fp for fp in folder_prefixes], second_basic_info_paths)):
        prefix = folder_prefixes[g_idx]  # used for output file naming
        log_file = out_dir / f"{folder_prefixes[g_idx]}_annotation.log"
        final_basic = run_mirn_rename(args, seconds, out_dir, prefix, log_file, project_root)
        basic_info_paths.append(final_basic)

    # ---- Consistency / common analysis (4.11) ----
    consistency_ref = args.consistency
    common_flag = args.common
    if consistency_ref or common_flag:
        # Use the final basic-info files (output/{prefix}-basic-info)
        info_files = basic_info_paths
        if num_groups == 1 and common_flag:
            print("Warning: --common requires multiple groups, skipping consistency step.")
        else:
            # ... same consistency logic as before but using info_files ...
            if consistency_ref and not common_flag:
                for info_path in info_files:
                    out_cons = info_path.parent / f"{info_path.stem}-consistency"
                    subprocess.run([sys.executable, str(project_root / "src" / "miRNA_consistency.py"),
                                    "-i", str(info_path), "-b", consistency_ref, "-ob", str(out_cons)],
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
                basic_info_paths = [p.parent / f"{p.stem}-consistency" for p in info_files]
            elif common_flag and not consistency_ref:
                inputs_str = ",".join(str(p) for p in info_files)
                outputs_common = [p.parent / f"{p.stem}-common" for p in info_files]
                out_str = ",".join(str(o) for o in outputs_common)
                subprocess.run([sys.executable, str(project_root / "src" / "miRNA_consistency.py"),
                                "-i", inputs_str, "--common", "--output", out_str],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
                basic_info_paths = outputs_common
            else:  # both
                inputs_str = ",".join(str(p) for p in info_files)
                outputs_cons = [p.parent / f"{p.stem}-consistency" for p in info_files]
                outputs_comm = [p.parent / f"{p.stem}-common" for p in info_files]
                cons_str = ",".join(str(o) for o in outputs_cons)
                comm_str = ",".join(str(o) for o in outputs_comm)
                subprocess.run([sys.executable, str(project_root / "src" / "miRNA_consistency.py"),
                                "-i", inputs_str, "-b", consistency_ref,
                                "-ob", cons_str, "--common", "--output", comm_str],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
                basic_info_paths = outputs_comm

    # ---- Execute part 2 (cluster, expression) for each group ----
    for g_idx, grp in enumerate(replicate_groups):
        out_dir = output_root / folder_prefixes[g_idx]
        group_fprefixes = [file_prefixes[input_folders.index(inp)] for inp in grp]
        print(f"  [finish] group {g_idx+1} ({folder_prefixes[g_idx]}) finalizing...")
        run_annotation_part2(args, grp, group_fprefixes, out_dir, genome, basic_info_paths[g_idx], project_root)

    # ---- Clean temporary directories if requested ----
    if args.clean:
        print("Cleaning temporary directories...")
        for out_dir in [output_root / fp for fp in folder_prefixes]:
            for sub in ['temp', 'index']:
                target = out_dir / sub
                if target.exists():
                    shutil.rmtree(target)

    print("\nAnnotation step completed successfully.")
    
    
def build_parser():
    """Return an independent parser for annotation arguments (without subcommand)."""
    parser = argparse.ArgumentParser(add_help=False)
    add_arguments(parser)
    return parser
