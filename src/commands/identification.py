#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
MirDeep-P3 identification step orchestrator.
Handles input parsing, dependency checks, output directory structure,
pipeline parallelisation, and log file generation.
"""

import argparse
import os
import sys
import shutil
import subprocess
import multiprocessing
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional

from utils.config import load_config
from utils.dependencies import check_external_tools
from check_python_deps import check_python_deps
from datetime import datetime

# ----------------------------------------------------------------------
# Constants
DEFAULT_MIN_LEN = 18
DEFAULT_MAX_LEN = 26
DEFAULT_RPM_THRESHOLD = 5
DEFAULT_MAX_MAPPINGS = 15
DEFAULT_PRE_LENGTH = 300
DEFAULT_THREADS = 1
DEFAULT_PROGRESS = 1

# ----------------------------------------------------------------------
def add_arguments(parser: argparse.ArgumentParser):
    """Define all command-line arguments for the identification step."""
    # Input / output
    parser.add_argument("-i", "--input", help="Input FASTQ/FASTA file(s), comma separated")
    parser.add_argument("-f", "--file", help="File containing list of input files, one per line")
    parser.add_argument("-o", "--output", help="Output root directory")
    parser.add_argument("--prefix", help="Output prefix(es) comma separated")
    parser.add_argument("--type", choices=["fastq", "fq", "fasta", "fa"], help="Input file type")
    parser.add_argument("-g", "--genome", help="Reference genome FASTA file")
    parser.add_argument("-d", "--index", help="Bowtie index prefix (if provided, skip building)")

    # Replicate handling
    parser.add_argument("-r", "--replicate",
                        help="Replicate grouping: comma-separated counts per group "
                             "(default: all inputs as one group, i.e., -r equals number of input files)")

    # Resource
    parser.add_argument("-t", "--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("-p", "--progress", type=int, default=DEFAULT_PROGRESS,
                        help="Number of parallel processes")

    # Read processing
    parser.add_argument("--min-len", type=int, default=DEFAULT_MIN_LEN)
    parser.add_argument("--max-len", type=int, default=DEFAULT_MAX_LEN)
    parser.add_argument("--rpm-threshold", type=float, default=DEFAULT_RPM_THRESHOLD)
    parser.add_argument("--max-mappings", type=int, default=DEFAULT_MAX_MAPPINGS)
    parser.add_argument("--pre-length", type=int, default=DEFAULT_PRE_LENGTH)

    # External tool paths
    parser.add_argument("--trim_galore", help="path to trim_galore")
    parser.add_argument("--bowtie", help="path to bowtie")
    parser.add_argument("--bowtie-build", help="path to bowtie-build")
    parser.add_argument("--RNAfold", help="path to RNAfold")
    parser.add_argument("--bedtools", help="path to bedtools")
    parser.add_argument("--samtools", help="path to samtools")
    parser.add_argument("--clean", action="store_true",
                        help="Remove temporary directories (temp, index) after successful completion")

    # Switches
    parser.add_argument("--reads_clean", dest="reads_clean", action="store_true", default=True,
                        help="Enable read cleaning (default on)")
    parser.add_argument("--no-reads_clean", dest="reads_clean", action="store_false",
                        help="Disable read cleaning")

    parser.add_argument("-h", "--help", action="help", help="Show this help message and exit.")

# ----------------------------------------------------------------------
def parse_input_files(config: Dict, args) -> List[Path]:
    """
    Resolve input file list.
    Priority:
      1) --file / -f (read lines from file)
      2) --input / -i (comma-separated list)
      3) config i/input or input (file or comma-separated)
    Exits on error if no input found.
    """
    files_raw = []

    if args.file:
        # Read from file (one path per line, ignore empty/whitespace lines)
        file_path = Path(args.file)
        if not file_path.is_file():
            sys.exit(f"Input list file not found: {args.file}")
        with open(file_path) as f:
            files_raw = [line.strip() for line in f if line.strip()]
        if args.input:
            print("[info] Both --file and --input given; using --file and ignoring --input.")
    elif args.input:
        files_raw = [p.strip() for p in args.input.split(",")]
    else:
        # Try config
        raw_cfg = config.get("i/input") or config.get("input")
        if not raw_cfg:
            sys.exit("Error: input files must be specified via -i, -f, or config.")
        if Path(raw_cfg).is_file():
            with open(raw_cfg) as f:
                files_raw = [line.strip() for line in f if line.strip()]
        else:
            files_raw = [p.strip() for p in raw_cfg.split(",")]

    # Convert to Path and verify existence
    input_files = [Path(p) for p in files_raw if p]
    if not input_files:
        sys.exit("Error: no valid input file paths found.")
    for fp in input_files:
        if not fp.is_file():
            sys.exit(f"Input file not found: {fp}")
    return input_files

def auto_detect_type(file_path: Path) -> str:
    """Return 'fastq' or 'fasta' based on extension."""
    name = file_path.name.lower()
    if any(name.endswith(ext) for ext in ['.fastq', '.fq', '.fastq.gz', '.fq.gz']):
        return 'fastq'
    elif any(name.endswith(ext) for ext in ['.fasta', '.fa', '.fasta.gz', '.fa.gz']):
        return 'fasta'
    else:
        return 'fastq'  # default

def decompress_file(input_path: Path, temp_dir: Path) -> Path:
    """If input is .gz, decompress to temp_dir and return path to decompressed file."""
    if input_path.suffix == '.gz':
        out_path = temp_dir / input_path.stem  # removes .gz
        if not out_path.exists():
            subprocess.run(f"gzip -dc {input_path} > {out_path}", shell=True, check=True)
        return out_path
    return input_path

def resolve_prefixes(args, config, num_inputs, replicate_groups):
    """Determine output directory prefixes for each input file."""
    prefix_raw = args.prefix or config.get("prefix", "")
    if not prefix_raw:
        # Auto‑prefix: use filename without extension for each input
        def stem(p):
            name = p.name
            if name.endswith('.gz'):
                name = name[:-3]
            for ext in ('.fastq', '.fq', '.fasta', '.fa'):
                if name.endswith(ext):
                    name = name[:-len(ext)]
                    break
            return name
        prefixes = [stem(f) for f in args.input_files]
        if len(prefixes) != num_inputs:
            sys.exit("Auto‑prefix error: count mismatch")
        return prefixes

    # User provided prefix(es)
    prefixes = [p.strip() for p in prefix_raw.split(",")]

    # Validate total inputs vs group sizes
    group_sizes = [len(g) for g in replicate_groups]
    total_expected = sum(group_sizes)
    if total_expected != num_inputs:
        sys.exit(f"Error: sum of replicate counts ({total_expected}) does not match number of inputs ({num_inputs}).")

    if len(prefixes) == 1:
        # Expand one prefix per group
        expanded = []
        for g_idx, g_size in enumerate(group_sizes):
            base = prefixes[0] + (f"-{g_idx+1}" if len(group_sizes) > 1 else "")
            for n in range(1, g_size+1):
                expanded.append(f"{base}-{n}")
        prefixes = expanded
    elif len(prefixes) == len(replicate_groups):
        # One prefix per group
        expanded = []
        for g_idx, g_size in enumerate(group_sizes):
            for n in range(1, g_size+1):
                expanded.append(f"{prefixes[g_idx]}-{n}")
        prefixes = expanded
    elif len(prefixes) != num_inputs:
        sys.exit(f"Error: prefix count ({len(prefixes)}) must be 1, {len(replicate_groups)} (number of groups), or {num_inputs} (number of inputs).")
    return prefixes


def ensure_defaults(args):
    """Ensure all identification-specific attributes exist on args with safe defaults."""
    defaults = {
        'clean': False,
        'reads_clean': True,
        'threads': DEFAULT_THREADS,
        'progress': DEFAULT_PROGRESS,
        'min_len': DEFAULT_MIN_LEN,
        'max_len': DEFAULT_MAX_LEN,
        'rpm_threshold': DEFAULT_RPM_THRESHOLD,
        'max_mappings': DEFAULT_MAX_MAPPINGS,
        'pre_length': DEFAULT_PRE_LENGTH,
        'input': None,
        'output': None,
        'prefix': None,
        'type': None,
        'genome': None,
        'index': None,
        'replicate': None,
        'trim_galore': None,
        'bowtie': None,
        'bowtie_build': None,
        'RNAfold': None,
        'bedtools': None,
        'samtools': None,
        'file': None,
    }
    for attr, val in defaults.items():
        if not hasattr(args, attr):
            setattr(args, attr, val)
            
            
# ----------------------------------------------------------------------
def run_pipeline_for_input(args, input_file: Path, prefix: str, output_root: Path,
                           genome: Path, index_prefix: str, log_file: Path):
    temp_dir = output_root / prefix / "temp"
    index_dir = output_root / prefix / "index"
    log_dir = log_file.parent
    log_dir.mkdir(parents=True, exist_ok=True)

    out_prefix = output_root / prefix
    out_prefix.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(parents=True, exist_ok=True)
    index_dir.mkdir(parents=True, exist_ok=True)

    # local shortcuts
    null = subprocess.DEVNULL
    py_exe = sys.executable
    threads = args.threads
    min_len = args.min_len
    max_len = args.max_len
    rpm_threshold = args.rpm_threshold
    max_mappings = args.max_mappings
    pre_length = args.pre_length
    trim_galore = getattr(args, 'trim_galore', None) or 'trim_galore'
    bowtie = getattr(args, 'bowtie', None) or 'bowtie'
    bowtie_build = getattr(args, 'bowtie_build', None) or 'bowtie-build'
    RNAfold = getattr(args, 'RNAfold', None) or 'RNAfold'
    samtools = getattr(args, 'samtools', None) or 'samtools'
    bedtools = getattr(args, 'bedtools', None) or 'bedtools'

    # Decompress if needed
    work_file = decompress_file(input_file, temp_dir)
    seq_type = args.type or auto_detect_type(work_file)
    is_fasta = seq_type in ('fasta', 'fa')

    # ---- 4.4 Read cleaning & 4.5 Preprocess / Copy ----
    if is_fasta:
        processed_fa = out_prefix / f"{prefix}.processed.fa"
        shutil.copy2(work_file, processed_fa)
        with open(log_file, 'a') as lf:
            lf.write(f"[info] FASTA input: copied {work_file.name} to {processed_fa}\n")
    else:
        if args.reads_clean:
            cmd = (f"{trim_galore} --small_rna --length {min_len} --max_length {max_len} "
                   f"--dont_gzip --suppress_warn -j 8 -o {out_prefix} {work_file}")
            subprocess.run(cmd, shell=True, check=True,
                           stdout=null, stderr=open(log_file, 'a'))
            # predict trimmed output
            trimmed_base = work_file.stem
            while '.' in trimmed_base:
                trimmed_base = Path(trimmed_base).stem
            trimmed_file = out_prefix / f"{trimmed_base}_trimmed.fq"
            if not trimmed_file.exists():
                sys.exit(f"Trimmed file not found: {trimmed_file}")
        else:
            trimmed_file = work_file

        processed_fa = temp_dir / f"{Path(trimmed_file).stem}.processed.fa"
        preprocess_script = args.project_root / "src" / "preprocess_files.py"
        cmd = (f"{py_exe} {preprocess_script} -i {trimmed_file} -o {temp_dir}")
        subprocess.run(cmd, shell=True, check=True,
                       stdout=null, stderr=open(log_file, 'a'))
        if not processed_fa.exists():
            sys.exit(f"Processed FASTA not created: {processed_fa}")

    # 4.6 Alignment to ncRNA databases
    project_data = args.project_root / "data"
    rfam_index = project_data / "index" / "rfam_index"
    mature_index = project_data / "index" / "mature_index"
    rfam_aln = temp_dir / "rfam_reads.aln"
    mature_aln = temp_dir / "mature_reads.aln"
    subprocess.run(f"{bowtie} -v 0 {rfam_index} -f {processed_fa} > {rfam_aln} 2>> {log_file}",
                   shell=True, check=True)
    subprocess.run(f"{bowtie} -v 1 {mature_index} -f {processed_fa} > {mature_aln} 2>> {log_file}",
                   shell=True, check=True)

    # 4.7 Filter ncRNA reads
    all_fa = temp_dir / f"{prefix}.fa"
    filtered_fa = temp_dir / f"{prefix}-processed.fa"
    total_reads_file = out_prefix / f"{prefix}.total_reads"
    preproc_reads_script = args.project_root / "src" / "preprocess_reads.py"
    cmd = (f"{py_exe} {preproc_reads_script} {processed_fa} {rfam_aln} {mature_aln} "
           f"{rpm_threshold} {all_fa} {filtered_fa} {total_reads_file}")
    subprocess.run(cmd, shell=True, check=True,
                   stdout=null, stderr=open(log_file, 'a'))

    # 4.8 Map to genome (index building)
    if args.index:
        genome_index = args.index
    else:
        genome_index = index_dir / "genome_index"
        if not (genome_index.with_suffix('.1.ebwt').exists() or genome_index.with_suffix('.1.ebwtl').exists()):
            subprocess.run(f"{bowtie_build} -f {genome} {genome_index} >> {log_file} 2>&1",
                           shell=True, check=True)

    is_large = bool(list(Path(genome_index).parent.glob(Path(genome_index).name + "*.ebwtl")))
    large_flag = " --large-index" if is_large else ""
    processed_aln = temp_dir / f"{prefix}-processed.aln"
    subprocess.run(f"{bowtie} -a -v 0{large_flag} {genome_index} -f {filtered_fa} > {processed_aln} 2>> {log_file}",
                   shell=True, check=True)

    # 4.9 Convert to BST and filter
    convert_script = args.project_root / "src" / "convert_bowtie_to_blast.py"
    filter_script = args.project_root / "src" / "filter_alignments.py"
    bst_file = temp_dir / f"{prefix}-processed.bst"
    subprocess.run(f"{py_exe} {convert_script} {processed_aln} {all_fa} {genome} -o {bst_file} "
                   f">> {log_file} 2>&1", shell=True, check=True)
    filter_bst = temp_dir / f"{prefix}-processed-filter.bst"
    subprocess.run(f"{py_exe} {filter_script} {bst_file} -c {max_mappings} -o {filter_bst} "
                   f">> {log_file} 2>&1", shell=True, check=True)

    # 4.10 Excise precursors and fold
    excise_script = args.project_root / "src" / "excise_candidate.py"
    precursors_fa = temp_dir / f"{prefix}_precursors.fa"
    subprocess.run(f"{py_exe} {excise_script} {genome} {filter_bst} -l {pre_length} -o {precursors_fa} "
                   f">> {log_file} 2>&1", shell=True, check=True)
    precursors_struc = temp_dir / f"{prefix}_precursors.struc"
    subprocess.run(f"{RNAfold} --noPS -j{threads} {precursors_fa} > {precursors_struc} 2>> {log_file}",
                   shell=True, check=True)

    # 4.11 Extract reads with no ncRNA
    all_aln = temp_dir / f"{prefix}.aln"
    subprocess.run(f"{bowtie} -a -v 0{large_flag} {genome_index} -f {all_fa} > {all_aln} 2>> {log_file}",
                   shell=True, check=True)
    all_bst = temp_dir / f"{prefix}.bst"
    subprocess.run(f"{py_exe} {convert_script} {all_aln} {all_fa} {genome} -o {all_bst} >> {log_file} 2>&1",
                   shell=True, check=True)
    all_filter_bst = temp_dir / f"{prefix}-filter.bst"
    subprocess.run(f"{py_exe} {filter_script} {all_bst} -c {max_mappings} -o {all_filter_bst} >> {log_file} 2>&1",
                   shell=True, check=True)
    filtered_fa_final = temp_dir / f"{prefix}_filtered.fa"
    subprocess.run(f"{py_exe} {filter_script} {all_filter_bst} -b {all_fa} -o {filtered_fa_final} >> {log_file} 2>&1",
                   shell=True, check=True)

    # 4.12 Prepare reads signature file
    prec_index = index_dir / f"{prefix}_precursors"
    subprocess.run(f"{bowtie_build} -f {precursors_fa} {prec_index} >> {log_file} 2>&1",
                   shell=True, check=True)
    prec_aln = temp_dir / f"{prefix}_precursors.aln"
    subprocess.run(f"{bowtie} -a -v 0 {prec_index} -f {filtered_fa_final} > {prec_aln} 2>> {log_file}",
                   shell=True, check=True)
    prec_bst = temp_dir / f"{prefix}_precursors.bst"
    subprocess.run(f"{py_exe} {convert_script} {prec_aln} {filtered_fa_final} {precursors_fa} -o {prec_bst} >> {log_file} 2>&1",
                   shell=True, check=True)
    signatures_file = temp_dir / f"{prefix}_signatures"
    subprocess.run(f"sort +3 -25 {prec_bst} > {signatures_file} 2>> {log_file}",
                   shell=True, check=True)

    # 4.13 miRNA prediction
    mod_mirdp_script = args.project_root / "src" / "mod_miRDP.py"
    predictions_file = temp_dir / f"{prefix}_predictions"
    subprocess.run(f"{py_exe} {mod_mirdp_script} {signatures_file} {precursors_struc} -o {predictions_file} >> {log_file} 2>&1",
                   shell=True, check=True)

    # 4.14 Plant-specific filtering
    fai_file = temp_dir / f"{genome.name}.fai"
    chrom_length_file = out_prefix / "chr_length"
    subprocess.run(f"{samtools} faidx {genome} --fai-idx {fai_file} > /dev/null 2>> {log_file}",
                   shell=True, check=True)
    with open(log_file, 'a') as lf:
        lf.write(f"Extracting chromosome lengths...\n")
    with open(fai_file) as fin, open(chrom_length_file, 'w') as fout:
        for line in fin:
            parts = line.split('\t')
            if len(parts) >= 2:
                fout.write(f"{parts[0]}\t{parts[1]}\n")
    nr_pred = out_prefix / f"{prefix}_nr_predictions"
    filter_pred = out_prefix / f"{prefix}_filter_P_prediction"
    rm_script = args.project_root / "src" / "mod_rm_redundant_meet_plant.py"
    subprocess.run(f"{py_exe} {rm_script} {chrom_length_file} {precursors_fa} {predictions_file} {total_reads_file} "
                   f"-n {nr_pred} -f {filter_pred} >> {log_file} 2>&1",
                   shell=True, check=True)

    # 4.15 Convert to BED
    bed_script = args.project_root / "src" / "convert_to_bed.py"
    bed_file = out_prefix / f"{prefix}_nr_predictions.bed"
    subprocess.run(f"{py_exe} {bed_script} -i {nr_pred} -o {bed_file} >> {log_file} 2>&1",
                   shell=True, check=True)
    
    # Copy filtered processed.fa to output directory for downstream use
    final_processed = temp_dir / f"{prefix}-processed.fa"
    if final_processed.exists():
        dest = out_prefix / f"{prefix}-processed.fa"
        shutil.copy2(final_processed, dest)
        with open(log_file, 'a') as lf:
            lf.write(f"[info] Copied {final_processed} to {dest}\n")

    return 0

def process_worker(args, input_file, prefix, root, genome, index, logf):
    """Wrapper that prints start/end for a single input file."""
    print(f"  [start] {prefix}")
    run_pipeline_for_input(args, input_file, prefix, root, genome, index, logf)
    print(f"  [done]  {prefix}")
# ----------------------------------------------------------------------
def run(args):
    """Main entry point for identification subcommand."""
    # 1. Load config and ensure defaults
    ensure_defaults(args)
    cfg = getattr(args, 'config_data', {})

    # 2. Combine config values into args (command line takes precedence)
    config_mapping = {
        'i/input': 'input',
        'o/output': 'output',
        'g/genome': 'genome',
        't/threads': 'threads',
        'p/progress': 'progress',
        'r/replicate': 'replicate',
        'd/index': 'index',
        'prefix': 'prefix',
        'min-length': 'min_len',
        'max-length': 'max_len',
        'type': 'type',
        'min-len': 'min_len',
        'max-len': 'max_len',
        'rpm-threshold': 'rpm_threshold',
        'max-mappings': 'max_mappings',
        'pre-length': 'pre_length',
        'trim_galore': 'trim_galore',
        'bowtie': 'bowtie',
        'bowtie-build': 'bowtie_build',
        'RNAfold': 'RNAfold',
        'bedtools': 'bedtools',
        'samtools': 'samtools',
        'reads_clean': 'reads_clean',
        'clean': 'clean',
        'qc': 'reads_clean',
    }
    for cfg_key, attr in config_mapping.items():
        if cfg_key in cfg and getattr(args, attr, None) is None:
            setattr(args, attr, cfg[cfg_key])

    # 3. Dependency detection
    required_tools = {
        'bowtie': getattr(args, 'bowtie', None) or cfg.get('bowtie'),
        'bowtie-build': getattr(args, 'bowtie_build', None) or cfg.get('bowtie-build'),
        'RNAfold': getattr(args, 'RNAfold', None) or cfg.get('RNAfold'),
        'samtools': getattr(args, 'samtools', None) or cfg.get('samtools'),
    }
    if getattr(args, 'reads_clean', True):
        required_tools['trim_galore'] = getattr(args, 'trim_galore', None) or cfg.get('trim_galore')
    # Optional tools for filtering and processing
    for opt_tool in ['bedtools', 'seqkit', 'cutadapt']:
        val = getattr(args, opt_tool, None) or cfg.get(opt_tool)
        if val:
            required_tools[opt_tool] = val

    all_ok, missing = check_external_tools(required_tools)
    if not all_ok:
        sys.exit("Error: missing required external dependencies.")
    check_python_deps()

    # ---- 2. Parameter resolution ----
    # Input files
    input_files = parse_input_files(cfg, args)
    if not input_files:
        sys.exit("Error: no input files provided.")
    args.input_files = input_files

    # Output root directory
    output_root = args.output or cfg.get("o/output") or cfg.get("output")
    if not output_root:
        step_str = "identification"
        timestamp = datetime.now().strftime("%m%d%y-%H%M")
        output_root = f"mirdeep-{step_str}-{timestamp}"
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    # Genome
    genome = args.genome or cfg.get("g/genome") or cfg.get("genome")
    if not genome:
        sys.exit("Error: genome file must be specified.")
    genome = Path(genome)
    if not genome.exists():
        sys.exit(f"Genome file not found: {genome}")

    # Replicate grouping
    rep_raw = args.replicate or cfg.get("r/replicate") or cfg.get("replicate") or str(len(input_files))
    rep_counts = [int(x.strip()) for x in rep_raw.split(",")]
    # Build groups
    replicate_groups = []
    idx = 0
    for cnt in rep_counts:
        replicate_groups.append(input_files[idx:idx+cnt])
        idx += cnt
    if idx != len(input_files):
        sys.exit("Error: replicate counts do not sum to number of input files.")
    args.replicate_groups = replicate_groups

    # Prefixes
    prefixes = resolve_prefixes(args, cfg, len(input_files), replicate_groups)
    args.prefixes = prefixes

    # Additional parameters (with defaults)
    args.min_len = args.min_len or cfg.get("min-len", DEFAULT_MIN_LEN)
    args.max_len = args.max_len or cfg.get("max-len", DEFAULT_MAX_LEN)
    args.rpm_threshold = args.rpm_threshold or cfg.get("rpm-threshold", DEFAULT_RPM_THRESHOLD)
    args.max_mappings = args.max_mappings or cfg.get("max-mappings", DEFAULT_MAX_MAPPINGS)
    args.pre_length = args.pre_length or cfg.get("pre-length", DEFAULT_PRE_LENGTH)
    args.threads = args.threads or cfg.get("t/threads", DEFAULT_THREADS)
    args.progress = args.progress or cfg.get("p/progress", DEFAULT_PROGRESS)

    print(f"\nProcessing {len(input_files)} input file(s) with {args.progress} parallel process(es).")
    # ---- 4. Create pipe file (overview) ----
    # output_root/{first_prefix?}_identification.pipe   Actually use a single pipe file per whole step.
    # We'll write it in output root.
    timestamp = datetime.now().strftime("%m%d%Y-%H%M")
    pipe_file = output_root / f"mirdp3-identification-{timestamp}.pipe"
    with open(pipe_file, 'w') as pf:
        for inp, pref in zip(input_files, prefixes):
            out_dir = output_root / pref
            group_prefix = pref.rsplit('-', 1)[0] if '-' in pref else pref
            pf.write(f"{inp}\t{out_dir}\t{group_prefix}\n")

    # ---- 5. Run pipelines with status tracking ----
    log_files = []
    for pref in prefixes:
        log_path = output_root / pref / f"{pref}_identification.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_files.append(log_path)

    if args.progress > 1 and len(input_files) > 1:
        with multiprocessing.Pool(processes=min(args.progress, len(input_files))) as pool:
            results = []
            for inp, pref, logf in zip(input_files, prefixes, log_files):
                results.append(pool.apply_async(process_worker,
                                                (args, inp, pref, output_root, genome, args.index, logf)))
            for res in results:
                res.get()
    else:
        for inp, pref, logf in zip(input_files, prefixes, log_files):
            process_worker(args, inp, pref, output_root, genome, args.index, logf)

    # ---- Clean temp file ----
    if args.clean:
        print("Cleaning temporary directories...")
        for pref in set(prefixes):
            for subdir in ['temp', 'index']:
                target = output_root / pref / subdir
                if target.exists():
                    shutil.rmtree(target)

    print("\nIdentification step completed successfully.")
    
def build_parser():
    """Return an independent parser for identification arguments (without subcommand)."""
    parser = argparse.ArgumentParser(add_help=False)
    add_arguments(parser)
    return parser