#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
Target_finder: miRNA target prediction.
Supports input via -i (FASTA) or -b (basic-info file from annotation).

Two interchangeable engines (--engine):
  * psrna (default): built-in psRNATarget-compatible scanner
    (src/MiRNATarget/psrna_scan.py). No external aligner, exhaustively optimises
    every site, and reproduces the psRNATarget web output far more closely
    (99.5% of its miRNA-target pairs vs 62.6% for ssearch on Arctium_lappa).
  * ssearch: ssearch36 + the MiRNATarget parser (the historical engine). Kept for
    backwards compatibility and for reproducing pre-3.1.7b results; its recall is
    limited by the ssearch36 -b HSP cap and by its heuristic coordinate parsing.

Both engines write the same <output>/target_finder.tsv with six columns:
    miRNA  target  mature_length  target_start  target_end  expectation
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

TARGET_ENGINES = ("ssearch", "psrna")
DEFAULT_ENGINE = "psrna"


def add_arguments(parser: argparse.ArgumentParser):
    """Define Target_finder-specific arguments."""
    parser.add_argument("-i", "--input", help="Input mature miRNA FASTA file")
    parser.add_argument("-b", "--basic", help="Input basic-info file from annotation (instead of -i)")
    parser.add_argument("-c", "--cds", required=True, help="CDS sequences FASTA file")
    parser.add_argument("-o", "--output", help="Output directory / prefix (def. mirdeep-target_finder-<timestamp>)")
    parser.add_argument("-e", "--evalue", type=float, default=2.5,
                        help="E-value cutoff for target filtering (default: 2.5)")
    parser.add_argument("--GUs", type=int, default=1,
                        help="Maximum number of G:U wobbles (default: 1; ssearch engine only)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads (default: 1)")
    parser.add_argument("--engine", choices=list(TARGET_ENGINES), default=DEFAULT_ENGINE,
                        help="Alignment engine: psrna (default, built-in psRNATarget-compatible "
                             "scanner, no ssearch36 needed) or ssearch (ssearch36 + MiRNATarget "
                             "parser, the historical engine)")
    parser.add_argument("--allow-gap", action=argparse.BooleanOptionalAction, default=None,
                        help="[psrna only] report single-base-indel sites too (default: on; "
                             "use --no-allow-gap for the indel-free subset only)")
    parser.add_argument("--max-hsp", type=int, default=200,
                        help="[ssearch only] ssearch36 -b value: maximum HSPs reported per "
                             "miRNA (default: 200, the historical behaviour). Raise it when a "
                             "miRNA has many targets and hits are being truncated.")
    parser.add_argument("--psrna-detail", action="store_true",
                        help="[psrna only] also write the native psRNATarget 14-column table "
                             "(<output>/target_finder.psrna.tsv)")


def _clean_name(name: str) -> str:
    """Strip trailing punctuation, matching psrna_scan.clean_name()."""
    name = name.strip()
    while name and not name[-1].isalnum():
        name = name[:-1]
    return name


def _read_fasta_lengths(path: Path) -> dict:
    """Map sequence name -> length (U->T, multi-line sequences supported)."""
    lengths = {}
    cur = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                cur = _clean_name(line[1:].split()[0])
                lengths[cur] = 0
            elif cur is not None:
                lengths[cur] += len(line)
    return lengths


def _fasta_order(path: Path) -> dict:
    """Map sequence name -> 0-based position in the FASTA (used for stable output order)."""
    order = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = _clean_name(line[1:].split()[0])
                order[name] = len(order)
    return order


def _sort_tab6(tsv: Path, order: dict) -> None:
    """Rewrite a six-column table in a deterministic order (input order, target, start)."""
    rows = []
    for line in open(tsv):
        line = line.rstrip("\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        try:
            start = int(fields[3])
        except ValueError:
            continue
        rows.append((order.get(fields[0], len(order)), fields[1], start, line))
    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    with open(tsv, "w") as out:
        for row in rows:
            out.write(row[3] + "\n")


def _write_tab6_from_psrna(raw_tsv: Path, out_tsv: Path, lengths: dict, order: dict) -> tuple:
    """Convert the native psRNATarget table into the six-column Target_finder layout."""
    rows = []
    for line in open(raw_tsv):
        line = line.rstrip("\n")
        if not line or line.startswith("miRNA_Acc."):
            continue
        fields = line.split("\t")
        if len(fields) < 9:
            continue
        try:
            start = int(fields[6])
            end = int(fields[7])
        except ValueError:
            continue
        length = lengths.get(fields[0]) or len(fields[8].replace("-", ""))
        rows.append((fields[0], fields[1], int(length), start, end, fields[2]))
    rows.sort(key=lambda r: (order.get(r[0], len(order)), r[1], r[3]))
    with open(out_tsv, "w") as out:
        for row in rows:
            out.write("%s\t%s\t%d\t%d\t%d\t%s\n" % row)
    return len(rows), len({(r[0], r[1]) for r in rows})


def _run_ssearch_engine(args, mature_fa: Path, cds_fa: Path, output_dir: Path, temp_dir: Path,
                        project_root: Path) -> Path:
    """Original ssearch36 + MiRNATarget pipeline."""
    mi_rna_target_dir = project_root / "src" / "MiRNATarget"
    scripts_dir = project_root / "scripts"
    threads = args.threads

    # ---- Step 1: Reverse complement ----
    rp_fa = temp_dir / "mature.rp.fa"
    subprocess.run(f"seqkit seq {shq(mature_fa)} -r -p > {shq(rp_fa)}",
                   shell=True, check=True)

    # ssearch36 base options
    max_hsp = getattr(args, "max_hsp", 200) or 200
    ssearch_opts = (f"-f -8 -g -3 -E 10000 -T {threads} -b {max_hsp} "
                    "-r +4/-3 -n -U -W 10 -N 20000")

    # Forward alignment
    forw_aln = temp_dir / "forward_alignment"
    subprocess.run(f"ssearch36 {ssearch_opts} {shq(mature_fa)} {shq(cds_fa)} "
                   f"> {shq(forw_aln)}",
                   shell=True, check=True)

    # Reverse alignment
    rev_aln = temp_dir / "reverse_alignment"
    subprocess.run(f"ssearch36 {ssearch_opts} {shq(rp_fa)} {shq(cds_fa)} "
                   f"> {shq(rev_aln)}",
                   shell=True, check=True)

    # ---- Step 3: Parse alignments to TSV ----
    parse_script = scripts_dir / "parse_ssearch.py"
    forw_tsv = temp_dir / "forward_alignment.tsv"
    rev_tsv = temp_dir / "reverse_alignment.tsv"
    subprocess.run(f"python {shq(parse_script)} -i {shq(forw_aln)} > {shq(forw_tsv)}",
                   shell=True, check=True)
    subprocess.run(f"python {shq(parse_script)} -i {shq(rev_aln)} > {shq(rev_tsv)}",
                   shell=True, check=True)
    total_tsv = temp_dir / "total.tsv"
    subprocess.run(f"cat {shq(forw_tsv)} {shq(rev_tsv)} > {shq(total_tsv)}",
                   shell=True, check=True)

    # ---- Step 4: Statistical filtering ----
    parse_targets_script = mi_rna_target_dir / "parse_mirna_targets.py"
    raw_target = temp_dir / "target_finder_raw.tsv"
    cmd = (f"python {shq(parse_targets_script)} -i {shq(total_tsv)} "
           f"--E_cutoff {args.evalue} --GUs_cutoff {args.GUs} > {shq(raw_target)}")
    subprocess.run(cmd, shell=True, check=True)

    # Get mature length
    mature_len = temp_dir / "mature.length"
    subprocess.run(f"seqkit fx2tab --length --name --header-line {shq(mature_fa)} | "
                   f"grep -v '#' > {shq(mature_len)}",
                   shell=True, check=True)

    # Final filtered output
    final_output = output_dir / "target_finder.tsv"
    subprocess.run(
        f"cat {shq(raw_target)} | grep -v '#' | "
        f"awk 'NR==FNR {{a[$1]=$2;next}} {{if($1 in a){{ print $0\"\\t\"a[$1] }}}}' "
        f"{shq(mature_len)} - | "
        f"awk '{{if($(NF-6)>=($NF-1)){{ print $0 }}}}' | "
        f"awk -F \"\\t\" '{{ print $1\"\\t\"$2\"\\t\"$NF\"\\t\"$8\"\\t\"($8+$NF-1)\"\\t\"$(NF-2) }}' "
        f"> {shq(final_output)}",
        shell=True, check=True)
    _sort_tab6(final_output, _fasta_order(mature_fa))
    return final_output


def _run_psrna_engine(args, mature_fa: Path, cds_fa: Path, output_dir: Path, temp_dir: Path,
                      project_root: Path) -> Path:
    """Built-in psRNATarget-compatible scanner."""
    scanner = project_root / "src" / "MiRNATarget" / "psrna_scan.py"
    if not scanner.is_file():
        sys.exit(f"Error: psrna engine not found: {scanner}")
    if args.GUs is not None and args.GUs != 1:
        print("[psrna] note: --GUs is ignored by the psrna engine "
              "(psRNATarget V2 has no G:U cap)")

    raw_tsv = temp_dir / "target_finder_raw.psrna.tsv"
    cmd = [sys.executable, str(scanner),
           "-m", str(mature_fa), "-t", str(cds_fa), "-o", str(raw_tsv),
           "-e", str(args.evalue), "-p", str(max(1, args.threads))]
    if args.allow_gap is not False:
        cmd.append("--allow-gap")
    subprocess.run(cmd, check=True)

    final_output = output_dir / "target_finder.tsv"
    n_sites, n_pairs = _write_tab6_from_psrna(raw_tsv, final_output,
                                              _read_fasta_lengths(mature_fa),
                                              _fasta_order(mature_fa))
    print(f"[psrna] {n_pairs} miRNA-target pairs, {n_sites} sites")
    if args.psrna_detail:
        detail = output_dir / "target_finder.psrna.tsv"
        shutil.copyfile(raw_tsv, detail)
        print(f"[psrna] native table written to {detail}")
    return final_output


def run(args):
    """Execute Target_finder analysis."""
    # ---- Resolve project root ----
    project_root = getattr(args, 'project_root', None) or Path(__file__).resolve().parents[2]

    # ---- Determine input ----
    if args.input and args.basic:
        sys.exit("Error: please specify either -i (FASTA) or -b (basic-info), not both.")
    if not args.input and not args.basic:
        sys.exit("Error: you must provide either -i or -b.")

    engine = getattr(args, "engine", DEFAULT_ENGINE) or DEFAULT_ENGINE
    if engine not in TARGET_ENGINES:
        sys.exit(f"Error: unknown engine '{engine}' (choose from {', '.join(TARGET_ENGINES)}).")
    if engine == "psrna" and getattr(args, "max_hsp", 200) != 200:
        sys.exit("Error: --max-hsp only applies to --engine ssearch.")
    if engine == "ssearch" and (args.allow_gap is not None or args.psrna_detail):
        sys.exit("Error: --allow-gap / --no-allow-gap / --psrna-detail only apply to "
                 "--engine psrna.")

    # ---- Output directory ----
    if args.output:
        output_dir = Path(args.output)
    else:
        output_dir = Path(f"mirdeep-target_finder-{datetime.now().strftime('%m%d%y-%H%M')}")
    output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = output_dir / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)

    # ---- Prepare mature FASTA ----
    if args.basic:
        basic_info_file = Path(args.basic).resolve()
        if not basic_info_file.is_file():
            sys.exit(f"Basic-info file not found: {basic_info_file}")
        mature_fa = temp_dir / "mature.fa"
        # Extract mature sequences: column 1 (ID) and column 18 (sequence)
        cmd = (f"cat {shq(basic_info_file)} | "
               f"awk '{{ print \">\"$1\"\\n\"$18 }}' > {shq(mature_fa)}")
        subprocess.run(cmd, shell=True, check=True)
        print(f"Extracted mature sequences from {basic_info_file} into {mature_fa}")
    else:
        mature_fa = Path(args.input).resolve()
        if not mature_fa.is_file():
            sys.exit(f"Input FASTA file not found: {mature_fa}")

    cds_fa = Path(args.cds).resolve()
    if not cds_fa.exists():
        sys.exit(f"CDS file not found: {cds_fa}")

    # ---- Check required external tools ----
    if engine == "ssearch":
        tools_ok, missing = check_external_tools({'ssearch36': None, 'seqkit': None})
        if not tools_ok:
            sys.exit("Error: ssearch36 / seqkit not found in PATH "
                     "(or use --engine psrna, which needs no external aligner).")

    # ---- Run ----
    if engine == "psrna":
        final_output = _run_psrna_engine(args, mature_fa, cds_fa, output_dir, temp_dir, project_root)
    else:
        final_output = _run_ssearch_engine(args, mature_fa, cds_fa, output_dir, temp_dir, project_root)

    print(f"\nTarget_finder completed. Results written to {final_output}")
