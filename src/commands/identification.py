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
import gzip
import os
import sys
import shutil
import subprocess
import multiprocessing
import zlib
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional
from collections import Counter

from utils.config import load_config
from utils.dependencies import check_external_tools
from utils.shellquote import shq
from utils.data_index import describe_index, index_status, installed_release
from check_python_deps import check_python_deps
from commands.stage_runner import (
    RunManifest,
    Stage,
    StageEmpty,
    digest,
    file_fingerprint,
    md5_file,
    run_stages,
)
from version import __version__ as MIRDP3_VERSION
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

# Per-sample outcome.  A sample that simply finds no candidate is a legitimate
# (empty) result, not an error: it must neither abort the run nor take down the
# other samples processed in the same batch.
PIPELINE_OK = "ok"
PIPELINE_EMPTY = "empty"
PIPELINE_FAILED = "failed"

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Input pre-flight validation
#
# These checks are deliberately cheap: a FASTQ is only inspected up to
# INPUT_VALIDATION_MAX_RECORDS records, so even a multi-gigabyte file costs
# milliseconds.  They exist so that a broken input is rejected *before* the
# expensive steps (genome index build, read cleaning, alignment) start, and so
# that the failure reads as one clear sentence instead of a traceback from
# whichever downstream tool happened to trip over it.
INPUT_VALIDATION_MAX_RECORDS = 4000
INPUT_VALIDATION_MAX_LINES = 200000
INPUT_VALIDATION_FOOTER_BYTES = 65536

# Nucleotide alphabet (DNA/RNA, gaps and wildcards tolerated).
_NUC_CHARS = set("ACGTUNacgtun-.*")


class InputValidationError(Exception):
    """A broken input file, described in one user-facing sentence."""


def _is_gzip_path(path: Path) -> bool:
    """True when the file is gzip-compressed, decided by magic bytes."""
    try:
        with open(path, 'rb') as fh:
            return fh.read(2) == b'\x1f\x8b'
    except OSError:
        return False


def _open_input_text(path: Path):
    """Open a plain or gzip-compressed input file as text."""
    if _is_gzip_path(path):
        return gzip.open(path, 'rt', encoding='utf-8', errors='replace', newline='')
    return open(path, 'rt', encoding='utf-8', errors='replace', newline='')


def _scan_fastq(fh, path: Path, max_records: int) -> Dict:
    """Validate the 4-line FASTQ structure over the first `max_records` records."""
    lineno = 0
    rec = 0
    eof = False
    seq_bases = 0
    seq_odd = 0

    while rec < max_records:
        # Skip blank separator lines; a trailing blank line is not an error.
        while True:
            line = fh.readline()
            if line == '':
                eof = True
                break
            lineno += 1
            if line.strip():
                break
        if eof:
            break

        block = [line.rstrip('\r\n')]
        while len(block) < 4:
            line = fh.readline()
            if line == '':
                break
            lineno += 1
            block.append(line.rstrip('\r\n'))

        rec += 1
        first = lineno - len(block) + 1
        if len(block) < 4:
            raise InputValidationError(
                "truncated FASTQ record #%d: the file ends after %d of 4 lines "
                "(line %d)" % (rec, len(block), first))

        header, seq, plus, qual = block
        if not header.startswith('@'):
            raise InputValidationError(
                "invalid FASTQ record #%d (line %d): header must start with '@', "
                "found %r" % (rec, first, header[:40]))
        if not plus.startswith('+'):
            raise InputValidationError(
                "invalid FASTQ record #%d (line %d): the third line must start "
                "with '+', found %r" % (rec, first + 2, plus[:40]))
        if not seq:
            raise InputValidationError(
                "invalid FASTQ record #%d (line %d): empty sequence" % (rec, first))
        if len(seq) != len(qual):
            raise InputValidationError(
                "invalid FASTQ record #%d (line %d): sequence length (%d) does not "
                "match quality length (%d)" % (rec, first, len(seq), len(qual)))

        seq_bases += len(seq)
        seq_odd += sum(1 for ch in seq if ch not in _NUC_CHARS)

    if not rec:
        raise InputValidationError("no FASTQ records found")
    if seq_odd > 0.5 * seq_bases:
        raise InputValidationError(
            "sequence does not look like DNA/RNA: %d of %d characters in the first "
            "%d record(s) are not nucleotide characters" % (seq_odd, seq_bases, rec))
    return {'records': rec, 'eof': eof}


def _scan_fasta(fh, path: Path, max_records: int, max_lines: int) -> Dict:
    """Validate that a FASTA file really contains nucleotide sequence records."""
    lineno = 0
    rec = 0
    bases = 0
    odd = 0
    seq_lines = 0
    hdr_line = 0
    eof = False

    for raw in fh:
        lineno += 1
        if lineno > max_lines:
            break
        text = raw.rstrip('\r\n')
        if not text.strip():
            continue
        if text.startswith('>'):
            if rec and not seq_lines:
                raise InputValidationError(
                    "FASTA record #%d (header on line %d) has no sequence"
                    % (rec, hdr_line))
            rec += 1
            if rec > max_records:
                rec = max_records
                break
            hdr_line = lineno
            seq_lines = 0
            continue
        if not rec:
            raise InputValidationError(
                "not a FASTA file: the first non-blank line (line %d) does not start "
                "with '>', found %r" % (lineno, text[:40]))
        seq_lines += 1
        bases += len(text)
        odd += sum(1 for ch in text if ch not in _NUC_CHARS)
    else:
        eof = True

    if rec and not seq_lines:
        raise InputValidationError(
            "FASTA record #%d (header on line %d) has no sequence" % (rec, hdr_line))
    if not rec:
        raise InputValidationError("no FASTA records found (no '>' header in the file)")
    if not bases:
        raise InputValidationError("FASTA contains headers but no nucleotide sequence")
    if odd > 0.5 * bases:
        raise InputValidationError(
            "not a nucleotide FASTA: %d of %d characters are not nucleotide "
            "characters" % (odd, bases))
    return {'records': rec, 'eof': eof}


def _check_plain_fastq_footer(path: Path) -> None:
    """
    Verify that an uncompressed FASTQ file ends on a complete record.

    A file truncated by an interrupted copy keeps a perfectly valid head, so the
    structure scan alone cannot see it; the last four lines are the cheap tell.
    """
    size = path.stat().st_size
    with open(path, 'rb') as fh:
        fh.seek(max(0, size - INPUT_VALIDATION_FOOTER_BYTES))
        tail = fh.read().decode('utf-8', 'replace')
    lines = [l for l in tail.splitlines() if l.strip()]
    if len(lines) < 4:
        return                                  # short file: the scan covered it
    header, seq, plus, qual = lines[-4:]
    if not header.startswith('@') or not plus.startswith('+') or len(seq) != len(qual):
        raise InputValidationError(
            "the file does not end on a complete FASTQ record, so it looks "
            "truncated (last four lines: %r, %r, %r, %r)"
            % (header[:30], seq[:30], plus[:30], qual[:30]))


def _verify_gzip_stream(path: Path) -> None:
    """Read a whole gzip stream to EOF, discarding the output, to check integrity."""
    try:
        with gzip.open(path, 'rb') as fh:
            while fh.read(1 << 22):
                pass
    except (gzip.BadGzipFile, zlib.error, EOFError) as exc:
        raise InputValidationError(
            "corrupt or truncated gzip input (%s: %s)"
            % (type(exc).__name__, exc))


def validate_input_file(path: Path, seq_type: Optional[str] = None,
                        verify_gzip: bool = False) -> Dict:
    """
    Inspect one input file before any computation is started.

    Returns a short description dict on success; raises InputValidationError with
    a concise message otherwise.  A gzip stream is fully integrity-checked
    whenever the scan reaches EOF; with `verify_gzip` the whole stream is read so
    that even a tail-only truncation is caught here instead of later.
    """
    path = Path(path)
    if not path.exists():
        raise InputValidationError("input file not found: %s" % path)
    if path.is_dir():
        raise InputValidationError("input path is a directory, not a file: %s" % path)
    try:
        size = path.stat().st_size
    except OSError as exc:
        raise InputValidationError("cannot stat input file: %s (%s)"
                                   % (path, exc.strerror or exc))
    if size == 0:
        raise InputValidationError("input file is empty (0 bytes): %s" % path)

    kind = (seq_type or auto_detect_type(path)).lower()
    kind = 'fasta' if kind in ('fasta', 'fa') else 'fastq'
    gz = path.suffix.lower() == '.gz'
    if gz and not _is_gzip_path(path):
        raise InputValidationError(
            "file has a .gz extension but is not gzip-compressed (bad magic "
            "bytes): %s" % path)

    try:
        with _open_input_text(path) as fh:
            if kind == 'fasta':
                scan = _scan_fasta(fh, path, INPUT_VALIDATION_MAX_RECORDS,
                                   INPUT_VALIDATION_MAX_LINES)
            else:
                scan = _scan_fastq(fh, path, INPUT_VALIDATION_MAX_RECORDS)
        full_gzip = gz and (scan['eof'] or verify_gzip)
        if gz and verify_gzip and not scan['eof']:
            _verify_gzip_stream(path)
        if kind == 'fastq' and not gz and not scan['eof']:
            _check_plain_fastq_footer(path)
    except InputValidationError as exc:
        raise InputValidationError("%s -- %s" % (exc, path))
    except (gzip.BadGzipFile, zlib.error, EOFError) as exc:
        raise InputValidationError(
            "corrupt or truncated gzip input (%s: %s) -- %s"
            % (type(exc).__name__, exc, path))
    except OSError as exc:
        raise InputValidationError("cannot read input file (%s) -- %s"
                                   % (exc.strerror or exc, path))

    return {'seq_type': kind, 'records': scan['records'], 'eof': scan['eof'],
            'gzip': gz, 'gzip_full': full_gzip, 'size': size}


def describe_input(info: Dict) -> str:
    """One-line human summary of a successful validation."""
    kind = 'FASTA' if info['seq_type'] == 'fasta' else 'FASTQ'
    count = ("%d record(s)" % info['records'] if info['eof']
             else ">=%d record(s), scan window" % info['records'])
    parts = ["%s %.2f MB" % (kind, info['size'] / 1e6), count]
    if info['gzip']:
        parts.append("gzip integrity OK" if info.get('gzip_full')
                     else "gzip integrity verified on decompression")
    return ", ".join(parts)


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
    parser.add_argument("--verify-gzip", dest="verify_gzip", action="store_true",
                        help="Read every gzipped input to EOF during input validation, so a "
                             "corrupt or truncated .gz is rejected before any computation. "
                             "Off by default (the tail is then checked while decompressing).")
    parser.add_argument("--clean", action="store_true",
                        help="Remove the per-sample temporary directories (<prefix>/temp) after "
                             "successful completion. The genome index is kept; add --clean-index "
                             "to remove that as well.")
    parser.add_argument("--clean-index", dest="clean_index", action="store_true",
                        help="Also remove the auto-built genome index (<out>/index) when cleaning. "
                             "An index passed with -d is never touched.")

    parser.add_argument("--resume", action="store_true",
                        help="Resume an interrupted run: skip every stage whose declared "
                             "outputs already exist and whose recorded run signature still "
                             "matches (same input, genome, parameters and version). Implies "
                             "keeping temp/, so --clean is ignored when resuming.")

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
    """
    If input is .gz, decompress to temp_dir and return path to decompressed file.

    Decompression runs in-process so that a truncated or corrupt gzip stream is
    reported as a concise InputValidationError instead of a raw traceback from
    `gzip -dc`.  The output is written to a .part file and renamed only after the
    stream has been read to completion, so an interrupted earlier run can never
    leave a half-decompressed file behind that a later run silently reuses.
    """
    input_path = Path(input_path)
    if input_path.suffix.lower() != '.gz':
        return input_path

    out_path = temp_dir / input_path.stem          # removes .gz
    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path

    part_path = out_path.with_name(out_path.name + '.part')

    def _cleanup():
        for leftover in (part_path, out_path):
            try:
                leftover.unlink()
            except OSError:
                pass

    try:
        with gzip.open(input_path, 'rb') as fin, open(part_path, 'wb') as fout:
            shutil.copyfileobj(fin, fout, length=1 << 20)
    except (gzip.BadGzipFile, zlib.error, EOFError) as exc:
        _cleanup()
        raise InputValidationError(
            "corrupt or truncated gzip input (%s: %s) -- %s"
            % (type(exc).__name__, exc, input_path))
    except OSError as exc:
        _cleanup()
        raise InputValidationError(
            "cannot decompress input file (%s) -- %s"
            % (exc.strerror or exc, input_path))

    os.replace(part_path, out_path)
    return out_path

def _uniquify_prefixes(prefixes: List[str]) -> List[str]:
    """
    Make output prefixes unique while preserving input order.

    A name occurring more than once receives a '-1', '-2', ... suffix in order of
    appearance (leafAA, leafAA -> leafAA-1, leafAA-2).  Names that are already
    unique are left untouched.  A final pass guarantees a collision-free result
    even when the inputs already contain such suffixed names.
    """
    counts = Counter(prefixes)
    seen = Counter()
    out = []
    for p in prefixes:
        if counts[p] == 1:
            out.append(p)
        else:
            seen[p] += 1
            out.append(f"{p}-{seen[p]}")

    if len(set(out)) != len(out):
        final, used = [], set()
        for p in out:
            cand, k = p, 1
            while cand in used:
                k += 1
                cand = f"{p}-{k}"
            used.add(cand)
            final.append(cand)
        out = final

    if out != prefixes:
        dup = sorted({p for p in prefixes if counts[p] > 1})
        print(f"[info] duplicate input name(s) detected: {', '.join(dup)}")
        for a, b in zip(prefixes, out):
            if a != b:
                print(f"[info]   output prefix: {a} -> {b}")
    return out


def resolve_prefixes(args, config, num_inputs, replicate_groups):
    """
    Determine output directory prefixes for each input file.

    Returns:
        ``(prefixes, group_prefixes)`` -- two equal-length lists.
        ``group_prefixes[i]`` is the name of the replicate group that input *i*
        belongs to and is shared by every input of that group.  The group name is
        recorded while the prefixes are built, never recovered afterwards by
        splitting the prefix string (which would corrupt file names that
        legitimately contain a '-').
    """
    prefix_raw = args.prefix or config.get("prefix", "")

    if not prefix_raw:
        # Auto-prefix: use filename without extension for each input
        def stem(p):
            name = p.name
            if name.endswith('.gz'):
                name = name[:-3]
            for ext in ('.fastq', '.fq', '.fasta', '.fa'):
                if name.endswith(ext):
                    name = name[:-len(ext)]
                    break
            return name
        bases = [stem(f) for f in args.input_files]
        if len(bases) != num_inputs:
            sys.exit("Auto-prefix error: count mismatch")
        group_prefixes, idx = [], 0
        for grp in replicate_groups:
            group_prefixes.extend([bases[idx]] * len(grp))
            idx += len(grp)
        return _uniquify_prefixes(bases), group_prefixes

    # User provided prefix(es)
    prefixes = [p.strip() for p in prefix_raw.split(",")]

    # Validate total inputs vs group sizes
    group_sizes = [len(g) for g in replicate_groups]
    total_expected = sum(group_sizes)
    if total_expected != num_inputs:
        sys.exit(f"Error: sum of replicate counts ({total_expected}) does not match number of inputs ({num_inputs}).")

    if len(prefixes) == 1:
        # Expand one prefix per group
        expanded, group_prefixes = [], []
        for g_idx, g_size in enumerate(group_sizes):
            base = prefixes[0] + (f"-{g_idx+1}" if len(group_sizes) > 1 else "")
            for n in range(1, g_size+1):
                expanded.append(f"{base}-{n}")
                group_prefixes.append(base)
        prefixes = expanded
    elif len(prefixes) == len(replicate_groups):
        # One prefix per group
        expanded, group_prefixes = [], []
        for g_idx, g_size in enumerate(group_sizes):
            for n in range(1, g_size+1):
                expanded.append(f"{prefixes[g_idx]}-{n}")
                group_prefixes.append(prefixes[g_idx])
        prefixes = expanded
    elif len(prefixes) != num_inputs:
        sys.exit(f"Error: prefix count ({len(prefixes)}) must be 1, {len(replicate_groups)} (number of groups), or {num_inputs} (number of inputs).")
    else:
        # One prefix per input: the group is taken from its first input
        group_prefixes, idx = [], 0
        for grp in replicate_groups:
            group_prefixes.extend([prefixes[idx]] * len(grp))
            idx += len(grp)

    return _uniquify_prefixes(prefixes), group_prefixes

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
def is_bowtie_index_complete(index_prefix: Path) -> bool:
    """
    Check if all required bowtie index files exist.
    Supports both small (ebwt) and large (ebwtl) formats.
    """
    required_suffixes = ['.1', '.2', '.3', '.4', '.rev.1', '.rev.2']
    # Determine which extension is in use by checking for .1.ebwt or .1.ebwtl
    if Path(str(index_prefix) + '.1.ebwt').exists():
        ext = '.ebwt'
    elif Path(str(index_prefix) + '.1.ebwtl').exists():
        ext = '.ebwtl'
    else:
        return False

    for suffix in required_suffixes:
        if not Path(str(index_prefix) + suffix + ext).exists():
            return False
    return True


# ----------------------------------------------------------------------
def _is_empty_file(path: Path) -> bool:
    """True when the file is missing or has zero bytes."""
    try:
        return path.stat().st_size == 0
    except OSError:
        return True


def _write_empty_results(out_prefix: Path, prefix: str, log_file: Path, reason: str,
                         tag: str = "empty") -> None:
    """
    Emit the placeholder files a downstream step expects when a sample yields no
    candidate, so the sample still leaves behind a valid (empty) result.
    """
    with open(log_file, 'a') as lf:
        lf.write(f"[{tag}] {prefix}: {reason}\n")
    for name in (f"{prefix}_nr_predictions",
                 f"{prefix}_nr_predictions.bed",
                 f"{prefix}_filter_P_prediction"):
        target = out_prefix / name
        if not target.exists():
            target.touch()


def _log_line(log_file: Path, text: str) -> None:
    """Append one line to the pipeline log."""
    with open(log_file, 'a') as lf:
        lf.write(text if text.endswith("\n") else text + "\n")


class SampleContext:
    """Everything one sample's stages need: resolved parameters and paths."""

    def __init__(self, args, input_file: Path, prefix: str, output_root: Path,
                 genome: Path, genome_index: Path, log_file: Path):
        self.args = args
        self.input_file = Path(input_file)
        self.prefix = prefix
        self.output_root = Path(output_root)
        self.genome = Path(genome)
        self.genome_index = Path(genome_index)
        self.log_file = Path(log_file)
        self.project_root = Path(args.project_root)
        self.out_prefix = self.output_root / prefix
        self.temp_dir = self.out_prefix / "temp"

        self.threads = args.threads
        self.min_len = args.min_len
        self.max_len = args.max_len
        self.rpm_threshold = args.rpm_threshold
        self.max_mappings = args.max_mappings
        self.pre_length = args.pre_length
        self.trim_galore = getattr(args, 'trim_galore', None) or 'trim_galore'
        self.bowtie = getattr(args, 'bowtie', None) or 'bowtie'
        self.bowtie_build = getattr(args, 'bowtie_build', None) or 'bowtie-build'
        self.RNAfold = getattr(args, 'RNAfold', None) or 'RNAfold'
        self.samtools = getattr(args, 'samtools', None) or 'samtools'
        self.py_exe = sys.executable

        src = self.project_root / "src"
        self.preprocess_script = src / "preprocess_files.py"
        self.preprocess_reads_script = src / "preprocess_reads.py"
        self.convert_script = src / "convert_bowtie_to_blast.py"
        self.filter_script = src / "filter_alignments.py"
        self.excise_script = src / "excise_candidate.py"
        self.mod_mirdp_script = src / "mod_miRDP.py"
        self.rm_script = src / "mod_rm_redundant_meet_plant.py"
        self.bed_script = src / "convert_to_bed.py"

        self.work_file = self.input_file
        self.seq_type = getattr(args, 'type', None)
        self.is_fasta = False
        self.is_large_index = False

    # -- setup ---------------------------------------------------------
    def prepare(self) -> None:
        """Decompress the input (idempotent) and decide FASTA vs FASTQ."""
        self.out_prefix.mkdir(parents=True, exist_ok=True)
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        self.work_file = decompress_file(self.input_file, self.temp_dir)
        self.seq_type = self.args.type or auto_detect_type(self.work_file)
        self.is_fasta = self.seq_type in ('fasta', 'fa')
        self.is_large_index = bool(list(
            self.genome_index.parent.glob(self.genome_index.name + "*.ebwtl")))

    # -- artifact paths (declared by the stages, consulted by --resume) --
    def processed_reads(self) -> Path:
        """The read file the ncRNA alignment stage consumes."""
        if self.is_fasta:
            return self.out_prefix / (self.prefix + ".processed.fa")
        if self.args.reads_clean:
            return self.temp_dir / (self.prefix + "_trimmed.processed.fa")
        return self.temp_dir / (Path(self.work_file).stem + ".processed.fa")

    def filtered_reads(self) -> Path:
        """Reads left after the ncRNA screen (before the genome-wide filter)."""
        return self.temp_dir / (self.prefix + ".processed.fa")

    def temp_file(self, name: str) -> Path:
        return self.temp_dir / name

    def out_file(self, name: str) -> Path:
        return self.out_prefix / name

    def manifest_path(self) -> Path:
        return self.out_prefix / (self.prefix + ".run.json")

    def large_flag(self) -> str:
        """Literal bowtie flag; must stay unquoted so it expands into its own word."""
        return " --large-index" if self.is_large_index else ""


# ----------------------------------------------------------------------
# Stage bodies.  Each one performs exactly the work the previous straight-line
# implementation did for that section, and raises StageEmpty (instead of
# returning early) when the sample legitimately has no candidate.
# ----------------------------------------------------------------------
def _stage_preprocess(ctx: SampleContext) -> None:
    null = subprocess.DEVNULL
    if ctx.is_fasta:
        processed_fa = ctx.processed_reads()
        shutil.copy2(ctx.work_file, processed_fa)
        _log_line(ctx.log_file,
                  f"[info] FASTA input: copied {ctx.work_file.name} to {processed_fa}")
        return

    if ctx.args.reads_clean:
        # --basename pins the output name, so a complex input name cannot hide it
        cmd = (f"{shq(ctx.trim_galore)} --small_rna --length {ctx.min_len} "
               f"--max_length {ctx.max_len} "
               f"--dont_gzip --suppress_warn -j 8 --basename {shq(ctx.prefix)} "
               f"-o {shq(ctx.out_prefix)} {shq(ctx.work_file)}")
        subprocess.run(cmd, shell=True, check=True,
                       stdout=null, stderr=open(ctx.log_file, 'a'))
        trimmed_file = ctx.out_file(ctx.prefix + "_trimmed.fq")
        if not trimmed_file.exists():
            raise RuntimeError(f"trimmed file not found: {trimmed_file}")
        if _is_empty_file(trimmed_file):
            raise StageEmpty("no read survived read cleaning")
    else:
        trimmed_file = ctx.work_file

    processed_fa = ctx.processed_reads()
    cmd = (f"{shq(ctx.py_exe)} {shq(ctx.preprocess_script)} -i {shq(trimmed_file)} "
           f"-o {shq(ctx.temp_dir)}")
    subprocess.run(cmd, shell=True, check=True,
                   stdout=null, stderr=open(ctx.log_file, 'a'))
    if not processed_fa.exists():
        raise RuntimeError(f"processed FASTA not created: {processed_fa}")
    if _is_empty_file(processed_fa):
        raise StageEmpty("no read survived preprocessing")


def _stage_ncrna_alignment(ctx: SampleContext) -> None:
    project_data = ctx.project_root / "data"
    rfam_index = project_data / "index" / "rfam_index"
    mature_index = project_data / "index" / "mature_index"
    rfam_aln = ctx.temp_file("rfam_reads.aln")
    mature_aln = ctx.temp_file("mature_reads.aln")
    processed_fa = ctx.processed_reads()
    subprocess.run(f"{shq(ctx.bowtie)} -v 0 {shq(rfam_index)} -f {shq(processed_fa)} "
                   f"> {shq(rfam_aln)} 2>> {shq(ctx.log_file)}",
                   shell=True, check=True)
    subprocess.run(f"{shq(ctx.bowtie)} -v 1 {shq(mature_index)} -f {shq(processed_fa)} "
                   f"> {shq(mature_aln)} 2>> {shq(ctx.log_file)}",
                   shell=True, check=True)


def _stage_ncrna_filter(ctx: SampleContext) -> None:
    all_fa = ctx.temp_file(ctx.prefix + ".fa")
    filtered_fa = ctx.filtered_reads()
    total_reads_file = ctx.out_file(ctx.prefix + ".total_reads")
    cmd = (f"{shq(ctx.py_exe)} {shq(ctx.preprocess_reads_script)} "
           f"{shq(ctx.processed_reads())} "
           f"{shq(ctx.temp_file('rfam_reads.aln'))} "
           f"{shq(ctx.temp_file('mature_reads.aln'))} {ctx.rpm_threshold} "
           f"{shq(all_fa)} {shq(filtered_fa)} {shq(total_reads_file)}")
    subprocess.run(cmd, shell=True, check=True,
                   stdout=subprocess.DEVNULL, stderr=open(ctx.log_file, 'a'))


def _stage_genome_map(ctx: SampleContext) -> None:
    aln = ctx.temp_file(ctx.prefix + "-processed.aln")
    subprocess.run(f"{shq(ctx.bowtie)} -a -v 0{ctx.large_flag()} {shq(ctx.genome_index)} "
                   f"-f {shq(ctx.filtered_reads())} > {shq(aln)} "
                   f"2>> {shq(ctx.log_file)}",
                   shell=True, check=True)


def _stage_bst_filter(ctx: SampleContext) -> None:
    bst = ctx.temp_file(ctx.prefix + "-processed.bst")
    filter_bst = ctx.temp_file(ctx.prefix + "-processed-filter.bst")
    subprocess.run(f"{shq(ctx.py_exe)} {shq(ctx.convert_script)} "
                   f"{shq(ctx.temp_file(ctx.prefix + '-processed.aln'))} "
                   f"{shq(ctx.temp_file(ctx.prefix + '.fa'))} {shq(ctx.genome)} "
                   f"-o {shq(bst)} >> {shq(ctx.log_file)} 2>&1",
                   shell=True, check=True)
    subprocess.run(f"{shq(ctx.py_exe)} {shq(ctx.filter_script)} {shq(bst)} "
                   f"-c {ctx.max_mappings} -o {shq(filter_bst)} "
                   f">> {shq(ctx.log_file)} 2>&1", shell=True, check=True)
    if _is_empty_file(filter_bst):
        raise StageEmpty("no candidate reads after alignment filtering")


def _stage_excise_fold(ctx: SampleContext) -> None:
    precursors_fa = ctx.temp_file(ctx.prefix + "_precursors.fa")
    subprocess.run(f"{shq(ctx.py_exe)} {shq(ctx.excise_script)} {shq(ctx.genome)} "
                   f"{shq(ctx.temp_file(ctx.prefix + '-processed-filter.bst'))} "
                   f"-l {ctx.pre_length} -o {shq(precursors_fa)} "
                   f">> {shq(ctx.log_file)} 2>&1", shell=True, check=True)
    if _is_empty_file(precursors_fa):
        raise StageEmpty("no precursor sequence could be excised")
    precursors_struc = ctx.temp_file(ctx.prefix + "_precursors.struc")
    subprocess.run(f"{shq(ctx.RNAfold)} --noPS -j{ctx.threads} {shq(precursors_fa)} "
                   f"> {shq(precursors_struc)} 2>> {shq(ctx.log_file)}",
                   shell=True, check=True)


def _stage_genome_map_all(ctx: SampleContext) -> None:
    """Map *all* reads (ncRNA included), then keep the ones with no ncRNA hit."""
    all_fa = ctx.temp_file(ctx.prefix + ".fa")
    all_aln = ctx.temp_file(ctx.prefix + ".aln")
    all_bst = ctx.temp_file(ctx.prefix + ".bst")
    all_filter_bst = ctx.temp_file(ctx.prefix + "-filter.bst")
    filtered_fa_final = ctx.temp_file(ctx.prefix + "_filtered.fa")
    subprocess.run(f"{shq(ctx.bowtie)} -a -v 0{ctx.large_flag()} {shq(ctx.genome_index)} "
                   f"-f {shq(all_fa)} > {shq(all_aln)} 2>> {shq(ctx.log_file)}",
                   shell=True, check=True)
    subprocess.run(f"{shq(ctx.py_exe)} {shq(ctx.convert_script)} {shq(all_aln)} "
                   f"{shq(all_fa)} {shq(ctx.genome)} -o {shq(all_bst)} "
                   f">> {shq(ctx.log_file)} 2>&1", shell=True, check=True)
    subprocess.run(f"{shq(ctx.py_exe)} {shq(ctx.filter_script)} {shq(all_bst)} "
                   f"-c {ctx.max_mappings} -o {shq(all_filter_bst)} "
                   f">> {shq(ctx.log_file)} 2>&1", shell=True, check=True)
    subprocess.run(f"{shq(ctx.py_exe)} {shq(ctx.filter_script)} {shq(all_filter_bst)} "
                   f"-b {shq(all_fa)} -o {shq(filtered_fa_final)} "
                   f">> {shq(ctx.log_file)} 2>&1", shell=True, check=True)


def _stage_signatures(ctx: SampleContext) -> None:
    precursors_fa = ctx.temp_file(ctx.prefix + "_precursors.fa")
    prec_index = ctx.temp_file(ctx.prefix + "_precursors")
    prec_aln = ctx.temp_file(ctx.prefix + "_precursors.aln")
    prec_bst = ctx.temp_file(ctx.prefix + "_precursors.bst")
    signatures_file = ctx.temp_file(ctx.prefix + "_signatures")
    subprocess.run(f"{shq(ctx.bowtie_build)} -f {shq(precursors_fa)} {shq(prec_index)} "
                   f">> {shq(ctx.log_file)} 2>&1", shell=True, check=True)
    subprocess.run(f"{shq(ctx.bowtie)} -a -v 0 {shq(prec_index)} "
                   f"-f {shq(ctx.temp_file(ctx.prefix + '_filtered.fa'))} "
                   f"> {shq(prec_aln)} 2>> {shq(ctx.log_file)}",
                   shell=True, check=True)
    subprocess.run(f"{shq(ctx.py_exe)} {shq(ctx.convert_script)} {shq(prec_aln)} "
                   f"{shq(ctx.temp_file(ctx.prefix + '_filtered.fa'))} "
                   f"{shq(precursors_fa)} -o {shq(prec_bst)} "
                   f">> {shq(ctx.log_file)} 2>&1", shell=True, check=True)
    subprocess.run(f"sort +3 -25 {shq(prec_bst)} > {shq(signatures_file)} "
                   f"2>> {shq(ctx.log_file)}", shell=True, check=True)


def _stage_prediction(ctx: SampleContext) -> None:
    predictions_file = ctx.temp_file(ctx.prefix + "_predictions")
    subprocess.run(f"{shq(ctx.py_exe)} {shq(ctx.mod_mirdp_script)} "
                   f"{shq(ctx.temp_file(ctx.prefix + '_signatures'))} "
                   f"{shq(ctx.temp_file(ctx.prefix + '_precursors.struc'))} "
                   f"-o {shq(predictions_file)} >> {shq(ctx.log_file)} 2>&1",
                   shell=True, check=True)
    if _is_empty_file(predictions_file):
        raise StageEmpty("miRDP2 predicted no miRNA candidate")


def _stage_plant_filter(ctx: SampleContext) -> None:
    fai_file = ctx.temp_file(ctx.genome.name + ".fai")
    chrom_length_file = ctx.out_file("chr_length")
    subprocess.run(f"{shq(ctx.samtools)} faidx {shq(ctx.genome)} --fai-idx {shq(fai_file)} "
                   f"> /dev/null 2>> {shq(ctx.log_file)}",
                   shell=True, check=True)
    _log_line(ctx.log_file, "Extracting chromosome lengths...")
    with open(fai_file) as fin, open(chrom_length_file, 'w') as fout:
        for line in fin:
            parts = line.split('\t')
            if len(parts) >= 2:
                fout.write(f"{parts[0]}\t{parts[1]}\n")

    nr_pred = ctx.out_file(ctx.prefix + "_nr_predictions")
    filter_pred = ctx.out_file(ctx.prefix + "_filter_P_prediction")
    rm_cmd = (f"{shq(ctx.py_exe)} {shq(ctx.rm_script)} {shq(chrom_length_file)} "
              f"{shq(ctx.temp_file(ctx.prefix + '_precursors.fa'))} "
              f"{shq(ctx.temp_file(ctx.prefix + '_predictions'))} "
              f"{shq(ctx.out_file(ctx.prefix + '.total_reads'))} "
              f"-n {shq(nr_pred)} -f {shq(filter_pred)} >> {shq(ctx.log_file)} 2>&1")
    rm_rc = subprocess.run(rm_cmd, shell=True).returncode
    if rm_rc != 0:
        # An empty result is detected below; a non-zero exit here means the step
        # itself failed (unreadable input), so say so instead of raising an opaque
        # CalledProcessError.
        raise RuntimeError(
            f"plant-specific filtering failed (exit {rm_rc}); see {ctx.log_file}")
    if _is_empty_file(nr_pred):
        raise StageEmpty("no candidate passed the plant-specific filter")


def _stage_bed(ctx: SampleContext) -> None:
    nr_pred = ctx.out_file(ctx.prefix + "_nr_predictions")
    bed_file = ctx.out_file(ctx.prefix + "_nr_predictions.bed")
    subprocess.run(f"{shq(ctx.py_exe)} {shq(ctx.bed_script)} -i {shq(nr_pred)} "
                   f"-o {shq(bed_file)} >> {shq(ctx.log_file)} 2>&1",
                   shell=True, check=True)


def _stage_finalize(ctx: SampleContext) -> None:
    """Publish the filtered reads under the canonical name annotation looks for."""
    final_processed = ctx.filtered_reads()
    if not final_processed.exists():
        return
    dest = ctx.out_file(ctx.prefix + ".processed.fa")
    shutil.copy2(final_processed, dest)
    _log_line(ctx.log_file, f"[info] Copied {final_processed} to {dest}")
    # A compatibility symlink keeps the legacy "<prefix>-processed.fa" spelling
    # working for external callers; annotation accepts either name.
    legacy = ctx.out_file(ctx.prefix + "-processed.fa")
    if not legacy.exists():
        try:
            legacy.symlink_to(dest.name)
        except OSError as exc:
            _log_line(ctx.log_file,
                      f"[warn] could not create compatibility symlink {legacy}: {exc}")


# ----------------------------------------------------------------------
# Stage table.  Order matters; artifacts drive both --resume and the manifest.
# ----------------------------------------------------------------------
def _preprocess_artifacts(ctx: SampleContext):
    artifacts = []
    if not ctx.is_fasta and ctx.args.reads_clean:
        artifacts.append(ctx.out_file(ctx.prefix + "_trimmed.fq"))
    artifacts.append(ctx.processed_reads())
    return artifacts


def _bowtie_index_artifacts(prefix: Path):
    """The six files bowtie-build leaves behind, for either index format."""
    suffixes = ('.1', '.2', '.3', '.4', '.rev.1', '.rev.2')
    for ext in ('.ebwt', '.ebwtl'):
        candidates = [Path(str(prefix) + s + ext) for s in suffixes]
        if all(p.exists() for p in candidates):
            return candidates
    return [Path(str(prefix) + s + '.ebwt') for s in suffixes]


def _signature_artifacts(ctx: SampleContext):
    return (_bowtie_index_artifacts(ctx.temp_file(ctx.prefix + "_precursors"))
            + [ctx.temp_file(ctx.prefix + "_precursors.aln"),
               ctx.temp_file(ctx.prefix + "_precursors.bst"),
               ctx.temp_file(ctx.prefix + "_signatures")])


def build_stages(ctx: SampleContext):
    """Return the ordered stage list for one sample."""
    return [
        Stage("4.4-4.5", "preprocess", _preprocess_artifacts, _stage_preprocess),
        Stage("4.6", "ncrna-alignment",
              lambda c: [c.temp_file("rfam_reads.aln"), c.temp_file("mature_reads.aln")],
              _stage_ncrna_alignment),
        Stage("4.7", "ncrna-filter",
              lambda c: [c.temp_file(c.prefix + ".fa"), c.filtered_reads(),
                         c.out_file(c.prefix + ".total_reads")],
              _stage_ncrna_filter),
        Stage("4.8", "genome-map",
              lambda c: [c.temp_file(c.prefix + "-processed.aln")],
              _stage_genome_map),
        Stage("4.9", "bst-filter",
              lambda c: [c.temp_file(c.prefix + "-processed.bst"),
                         c.temp_file(c.prefix + "-processed-filter.bst")],
              _stage_bst_filter),
        Stage("4.10", "excise-fold",
              lambda c: [c.temp_file(c.prefix + "_precursors.fa"),
                         c.temp_file(c.prefix + "_precursors.struc")],
              _stage_excise_fold),
        Stage("4.11", "genome-map-all",
              lambda c: [c.temp_file(c.prefix + ".aln"), c.temp_file(c.prefix + ".bst"),
                         c.temp_file(c.prefix + "-filter.bst"), c.filtered_reads()],
              _stage_genome_map_all),
        Stage("4.12", "signatures", _signature_artifacts, _stage_signatures),
        Stage("4.13", "prediction",
              lambda c: [c.temp_file(c.prefix + "_predictions")],
              _stage_prediction),
        Stage("4.14", "plant-filter",
              lambda c: [c.out_file("chr_length"),
                         c.out_file(c.prefix + "_nr_predictions"),
                         c.out_file(c.prefix + "_filter_P_prediction")],
              _stage_plant_filter),
        Stage("4.15", "bed",
              lambda c: [c.out_file(c.prefix + "_nr_predictions.bed")],
              _stage_bed),
        Stage("4.16", "finalize",
              lambda c: [c.out_file(c.prefix + ".processed.fa")],
              _stage_finalize),
    ]


def _manifest_signature(ctx: SampleContext) -> str:
    """Fingerprint of everything that would invalidate a resumed run."""
    return digest({
        "version": MIRDP3_VERSION,
        "prefix": ctx.prefix,
        # A data-index upgrade changes ncRNA masking, so it invalidates a resume (O-10).
        "data_index": installed_release(ctx.project_root / "data" / "index"),
        "input": file_fingerprint(ctx.input_file),
        "genome": file_fingerprint(ctx.genome),
        "genome_index": str(ctx.genome_index),
        "params": {
            "type": ctx.seq_type,
            "reads_clean": bool(ctx.args.reads_clean),
            "min_len": ctx.min_len,
            "max_len": ctx.max_len,
            "rpm_threshold": ctx.rpm_threshold,
            "max_mappings": ctx.max_mappings,
            "pre_length": ctx.pre_length,
        },
    })


def _manifest_params(ctx: SampleContext) -> Dict:
    return {
        "type": ctx.seq_type,
        "reads_clean": bool(ctx.args.reads_clean),
        "min_len": ctx.min_len,
        "max_len": ctx.max_len,
        "rpm_threshold": ctx.rpm_threshold,
        "max_mappings": ctx.max_mappings,
        "pre_length": ctx.pre_length,
        "threads": ctx.threads,
        "data_index": installed_release(ctx.project_root / "data" / "index"),
    }


def run_pipeline_for_input(args, input_file: Path, prefix: str, output_root: Path,
                           genome: Path, genome_index: Path, log_file: Path):
    ctx = SampleContext(args, input_file, prefix, output_root, genome, genome_index, log_file)
    ctx.prepare()

    stages = build_stages(ctx)
    resume = bool(getattr(args, 'resume', False))
    signature = _manifest_signature(ctx)
    manifest_file = ctx.manifest_path()

    manifest = None
    if resume and manifest_file.is_file():
        try:
            manifest = RunManifest.load(manifest_file)
        except (OSError, ValueError) as exc:
            _log_line(log_file, f"[resume] manifest unreadable ({exc}); starting from scratch")
            manifest = None
        if manifest is not None and not manifest.signature_matches(signature):
            _log_line(log_file, "[resume] input/parameters changed -- starting from scratch")
            print(f"  [resume] {prefix}: input or parameters changed, not resuming",
                  flush=True)
            manifest = None
    if manifest is None:
        manifest = RunManifest.create(
            manifest_file,
            tool="mirdeep-p3",
            version=MIRDP3_VERSION,
            step="identification",
            prefix=prefix,
            input_file=ctx.input_file,
            genome=ctx.genome,
            params=_manifest_params(ctx),
            tools=getattr(args, 'tool_versions', {}) or {},
            signature=signature,
            input_md5=md5_file(ctx.input_file),
        )
        _log_line(log_file, f"[manifest] {manifest_file}")

    try:
        run_stages(stages, ctx, manifest, resume=resume, log_file=log_file)
    except StageEmpty as exc:
        _write_empty_results(ctx.out_prefix, prefix, log_file, exc.reason)
        manifest.finish(PIPELINE_EMPTY, exc.reason)
        return PIPELINE_EMPTY, exc.reason
    except Exception as exc:
        manifest.finish(PIPELINE_FAILED, f"{type(exc).__name__}: {exc}")
        raise
    manifest.finish(PIPELINE_OK, "")
    return PIPELINE_OK, ""


def probe_tool_versions(args) -> Dict[str, str]:
    """Record the version line of every external tool the pipeline may call."""
    commands = {
        'bowtie': [getattr(args, 'bowtie', None) or 'bowtie', '--version'],
        'bowtie-build': [getattr(args, 'bowtie_build', None) or 'bowtie-build', '--version'],
        'RNAfold': [getattr(args, 'RNAfold', None) or 'RNAfold', '--version'],
        'samtools': [getattr(args, 'samtools', None) or 'samtools', '--version'],
    }
    if getattr(args, 'reads_clean', True):
        commands['trim_galore'] = [getattr(args, 'trim_galore', None) or 'trim_galore',
                                   '--version']
    versions = {}
    for name, cmd in commands.items():
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
        except (OSError, subprocess.SubprocessError):
            versions[name] = "unavailable"
            continue
        text = (result.stdout or result.stderr or "").strip()
        versions[name] = text.splitlines()[0].strip() if text else "unknown"
    return versions



def process_worker(args, input_file, prefix, root, genome, shared_genome_index, logf):
    """
    Run one sample and report its outcome.

    Exceptions are deliberately caught here: a single failing sample must not
    abort the whole batch -- in parallel mode a raised exception would propagate
    through AsyncResult.get() and kill the entire run.
    """
    import traceback
    print(f"  [start] {prefix}", flush=True)
    try:
        status, message = run_pipeline_for_input(args, input_file, prefix, root,
                                                 genome, shared_genome_index, logf)
    except Exception as exc:
        with open(logf, 'a') as lf:
            lf.write(f"\n[error] sample '{prefix}' failed: {exc}\n")
            lf.write(traceback.format_exc())
        print(f"  [fail]  {prefix}  ({type(exc).__name__}: {exc})", flush=True)
        return PIPELINE_FAILED, f"{type(exc).__name__}: {exc}"
    print(f"  [{status}]  {prefix}", flush=True)
    return status, message

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

    # Record the external tool versions once per run; every sample manifest
    # inherits them, so a result can be traced back to the binaries that made it.
    args.tool_versions = probe_tool_versions(args)
    if getattr(args, 'resume', False) and getattr(args, 'clean', False):
        # Resuming relies on the intermediate files kept in temp/.
        print("[resume] --clean ignored: temp/ is needed to resume.")
        args.clean = False

    # ---- 2. Parameter resolution ----
    # Input files
    input_files = parse_input_files(cfg, args)
    if not input_files:
        sys.exit("Error: no input files provided.")
    args.input_files = input_files

    # ---- 2c. Input pre-flight: reject broken files before any computation ----
    # A damaged file used to surface thirty seconds later as an opaque bowtie or
    # bowtie-build failure (or as a raw CalledProcessError traceback from gzip).
    # Everything that can be decided from the file itself is decided here, and a
    # rejected sample is skipped instead of taking down the whole batch.
    input_errors = {}
    print(f"\nValidating {len(input_files)} input file(s)...")
    for idx, fp in enumerate(input_files):
        try:
            info = validate_input_file(fp, getattr(args, 'type', None) or cfg.get('type'),
                                       verify_gzip=getattr(args, 'verify_gzip', False))
        except InputValidationError as exc:
            input_errors[idx] = str(exc)
            print(f"  [invalid] {exc}")
        else:
            print(f"  [ok]      {fp}  ({describe_input(info)})")

    if input_errors:
        if len(input_errors) == len(input_files):
            sys.exit("[error] every input file failed validation -- nothing to process.")
        print(f"[warn] {len(input_errors)} of {len(input_files)} input file(s) failed "
              f"validation and will be skipped; the remaining "
              f"{len(input_files) - len(input_errors)} will still be processed.")

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
    # Same pre-flight as the read files: a 0-byte or non-FASTA genome used to
    # abort with a raw CalledProcessError traceback from bowtie-build, which runs
    # outside the per-sample worker and therefore was not contained.
    try:
        ginfo = validate_input_file(genome, 'fasta')
    except InputValidationError as exc:
        sys.exit(f"[error] genome file is not usable: {exc}")
    print(f"[info] genome: {genome}  ({describe_input(ginfo)})")

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

    # Prefixes (group_prefixes records the replicate group of every input)
    prefixes, group_prefixes = resolve_prefixes(args, cfg, len(input_files), replicate_groups)
    args.prefixes = prefixes
    args.group_prefixes = group_prefixes

    # Additional parameters (with defaults)
    args.min_len = args.min_len or cfg.get("min-len", DEFAULT_MIN_LEN)
    args.max_len = args.max_len or cfg.get("max-len", DEFAULT_MAX_LEN)
    args.rpm_threshold = args.rpm_threshold or cfg.get("rpm-threshold", DEFAULT_RPM_THRESHOLD)
    args.max_mappings = args.max_mappings or cfg.get("max-mappings", DEFAULT_MAX_MAPPINGS)
    args.pre_length = args.pre_length or cfg.get("pre-length", DEFAULT_PRE_LENGTH)
    args.threads = args.threads or cfg.get("t/threads", DEFAULT_THREADS)
    args.progress = args.progress or cfg.get("p/progress", DEFAULT_PROGRESS)

    # Data index (O-10): software and index are released separately, so report
    # which index is installed instead of comparing it with the version.
    print(f"[index] data index: {describe_index(index_status(args.project_root))}")

    print(f"\nProcessing {len(input_files)} input file(s) with {args.progress} parallel process(es).")

    # ---- 4. Create pipe file (overview) ----
    timestamp = datetime.now().strftime("%m%d%Y-%H%M")
    # The batch manifest needs second resolution: two runs started in the same minute
    # must not overwrite each other's record (I-24).
    manifest_stamp = datetime.now().strftime("%m%d%Y-%H%M%S")
    pipe_file = output_root / f"mirdp3-identification-{timestamp}.pipe"
    with open(pipe_file, 'w') as pf:
        for inp, pref, gpref in zip(input_files, prefixes, group_prefixes):
            out_dir = output_root / pref
            pf.write(f"{inp}\t{out_dir}\t{gpref}\n")

    # ---- Build shared genome index (if needed) ----
    shared_index_dir = output_root / "index"
    shared_index_dir.mkdir(parents=True, exist_ok=True)
    if args.index:
        shared_genome_index = Path(args.index)
    else:
        shared_genome_index = shared_index_dir / "genome_index"
        if not is_bowtie_index_complete(shared_genome_index):
            print("Building index...")
            bowtie_build_exe = getattr(args, 'bowtie_build', None) or 'bowtie-build'
            index_log = output_root / "index_build.log"
            cmd = (f"{shq(bowtie_build_exe)} -f {shq(genome)} "
                   f"{shq(shared_genome_index)} >> {shq(index_log)} 2>&1")
            subprocess.run(cmd, shell=True, check=True)

    # ---- 5. Run pipelines with status tracking ----
    log_files = []
    for pref in prefixes:
        log_path = output_root / pref / f"{pref}_identification.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_files.append(log_path)

    outcomes = []          # (prefix, status, message)

    def _record(pref, ret):
        if isinstance(ret, tuple) and len(ret) == 2:
            outcomes.append((pref, ret[0], ret[1]))
        else:
            outcomes.append((pref, PIPELINE_OK, ""))

    # Samples rejected by the pre-flight pass are reported here and never
    # launched -- their inputs are known to be unreadable.
    runnable = []
    for idx, (inp, pref, logf) in enumerate(zip(input_files, prefixes, log_files)):
        if idx in input_errors:
            reason = f"invalid input: {input_errors[idx]}"
            _write_empty_results(output_root / pref, pref, logf, reason, tag="invalid")
            print(f"  [skip]  {pref}  (invalid input)")
            outcomes.append((pref, PIPELINE_FAILED, reason))
            continue
        runnable.append((inp, pref, logf))

    if args.progress > 1 and len(runnable) > 1:
        with multiprocessing.Pool(processes=min(args.progress, len(runnable))) as pool:
            jobs = []
            for inp, pref, logf in runnable:
                jobs.append((pref, pool.apply_async(
                    process_worker,
                    (args, inp, pref, output_root, genome, shared_genome_index, logf))))
            for pref, res in jobs:
                try:
                    _record(pref, res.get())
                except Exception as exc:      # defensive: worker already catches
                    outcomes.append((pref, PIPELINE_FAILED, f"{type(exc).__name__}: {exc}"))
    else:
        for inp, pref, logf in runnable:
            _record(pref, process_worker(args, inp, pref, output_root,
                                         genome, shared_genome_index, logf))

    # ---- 5b. Per-sample outcome summary (a bad sample never aborts the batch) ----
    ok_all = [p for p, s, _ in outcomes if s == PIPELINE_OK]
    empty  = [(p, m) for p, s, m in outcomes if s == PIPELINE_EMPTY]
    failed = [(p, m) for p, s, m in outcomes if s == PIPELINE_FAILED]

    print(f"\n[summary] {len(outcomes)} sample(s): {len(ok_all)} ok, "
          f"{len(empty)} empty (no candidate), {len(failed)} failed")
    if empty:
        print(f"[summary] EMPTY samples -- skipped, empty result written ({len(empty)}):")
        for p, m in empty:
            print(f"    - {p}  ({m})")
    if failed:
        print(f"[summary] FAILED samples ({len(failed)}):")
        for p, m in failed:
            print(f"    - {p}  ({m})")
        print("[summary] the remaining samples were still processed normally.")
    if failed and not ok_all and not empty:
        sys.exit("[error] every sample failed -- nothing was produced.")

    # ---- Clean temp file ----
    if args.clean:
        print("Cleaning temporary directories...")
        for pref in set(prefixes):
            for subdir in ['temp']:
                target = output_root / pref / subdir
                if target.exists():
                    shutil.rmtree(target)
    if getattr(args, 'clean_index', False):
        # Only the index this run built under <out>/index may be removed here; an
        # index supplied with -d belongs to the user and must be left alone.
        if getattr(args, 'index', None):
            print("[clean-index] skipped: the index came from -d and is left in place.")
        else:
            index_dir = output_root / "index"
            if index_dir.exists():
                print(f"[clean-index] removing {index_dir}")
                shutil.rmtree(index_dir)

    # ---- 5c. Batch manifest ----
    # One machine-readable record per batch, pointing at each sample manifest.
    batch_manifest = output_root / f"mirdp3-identification-{manifest_stamp}.manifest.json"
    batch = RunManifest.create_batch(
        batch_manifest,
        tool="mirdeep-p3",
        version=MIRDP3_VERSION,
        step="identification",
        params={
            "threads": args.threads,
            "progress": args.progress,
            "genome": str(genome),
            "genome_index": str(shared_genome_index),
            "resume": bool(getattr(args, 'resume', False)),
            "data_index": index_status(args.project_root),
        },
        tools=getattr(args, 'tool_versions', {}) or {},
        signature=digest({"prefixes": list(prefixes),
                          "genome": file_fingerprint(genome)}),
        inputs=input_files,
    )
    for pref, status, message in outcomes:
        batch.add_sample(pref, status, message,
                         str(output_root / pref / (pref + ".run.json")))
    if failed and (ok_all or empty):
        batch_status = "partial"
    elif failed:
        batch_status = PIPELINE_FAILED
    elif empty:
        batch_status = PIPELINE_EMPTY
    else:
        batch_status = PIPELINE_OK
    batch.finish(batch_status)
    print(f"[manifest] {batch_manifest}")

    print("\nIdentification step completed successfully.")

def build_parser():
    """Return an independent parser for identification arguments (without subcommand)."""
    parser = argparse.ArgumentParser(add_help=False)
    add_arguments(parser)
    return parser
