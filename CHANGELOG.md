# Changelog

All notable changes to MirDeep-P3 are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Release tags follow the scheme `v<version>` (source, used by the conda recipe) and
`v<version>-full` (Docker images + offline bundle, see `.github/workflows/docker-release.yml`).
Releases before 3.1.6c are listed on the
[GitHub Releases page](https://github.com/YangXZ-lab/mirdeep-p3/releases).

## [3.1.6c] - 2026-09-18

This release is about **not producing wrong or misleading results**. Three classes of
silent failure are fixed: damaged input that crashed a whole batch, legitimate empty
results that were reported as errors, and partial failures that were reported as success.

### Fixed

- **Minus-strand star coordinates** (`src/primary-basic-info.py`): the star sequence was
  placed at the wrong end of the precursor for every minus-strand record (171/171 in the
  test set had the star interval overlapping the mature miRNA). Coordinates are now
  mirrored according to `strand`, and the interval length matches the sequence length
  (was 20 nt for a 21 nt star). The duplicated logic in the conserved/non-conserved
  paths was replaced by a single `compute_star()`.

- **Empty result reported as an error** (`src/mod_rm_redundant_meet_plant.py`): a sample
  whose predictions could not be positioned (missing precursor, chromosome absent from
  the `chr_length` file, non-positive coordinates) exited 1, which aborted the run,
  produced a `CalledProcessError` at the top level and skipped the empty BED. Such a
  sample now writes empty outputs and exits 0, so it is recorded as an empty result.

- **False success in the enrichment chain** (`scripts/miRNA_enrich_analysis.py`,
  `src/commands/analysis/functional_analysis.py`, `src/commands/analysis/onestep.py`):
  per-miRNA failures were only printed, never aggregated, and the top level still
  printed "completed" and exited 0. Sub-task exit codes are now aggregated and every
  layer propagates the outcome (see *Changed* below).

- **`enrich_analysis.R` crashed on a NULL result**: `nrow(ekp)` and
  `lapply(ego_list, function(x) x@result)` aborted with
  `argument is of length zero` / a NULL-slot error when no gene could be mapped to a
  KEGG/GO term. A gene list with no significant term is now an empty result (exit code
  3), not a failure.

- **Damaged input was diagnosed far too late** (`src/commands/identification.py`):
  a truncated `.gz` produced a raw `CalledProcessError` traceback, while a truncated
  FASTQ or a 0-byte FASTA only failed ~30 s later inside bowtie — after the genome
  index had been built. Inputs are now checked in a pre-flight pass (gzip magic and
  integrity, FASTQ four-line structure, plain-FASTQ completeness, FASTA records,
  non-empty) before any computation starts, and the genome file is checked the same way.

- **Half-decompressed files were silently reused**: `decompress_file()` returned any
  existing output file, so a decompression interrupted by an earlier failure left a
  partial file that later runs consumed. Decompression now happens in-process, writes a
  `.part` file, and renames it atomically only after the stream has been read to the end.

- **`processed.fa` naming mismatch between steps**: identification copied the final
  filtered FASTA as `<prefix>-processed.fa` while annotation searched for
  `*.processed.fa`. FASTQ inputs therefore broke annotation. The canonical name is now
  `<prefix>.processed.fa` (with a compatibility symlink), and annotation accepts both
  spellings.

- **Paths containing spaces or parentheses broke nearly every command**: every path
  interpolated into a shell command string was unquoted, so the shell re-split it
  (e.g. `bowtie` exited 2 on `.../space test (v2)/leaf AA.fq`). All command
  constructions now quote their paths.

- **Duplicate input names collided**: two inputs named `leafAA.fastq` from different
  directories both wrote to `output/leafAA`. Duplicate names are now suffixed
  (`leafAA-1`, `leafAA-2`), and the replicate group name in the `.pipe` file is recorded
  when the prefix is built instead of being re-derived from the string.

- **Shared R library lock contention**: the OrgDb was installed into the shared conda R
  library, so concurrent jobs deadlocked on `00LOCK-*` and the base environment was
  polluted. Each run now installs into its own library (`MIRDEEP_R_LIB`), clears stale
  locks, and checks only that library before installing.

- **A single bad sample aborted the whole batch**: an exception in one sample propagated
  through `AsyncResult.get()` and killed the run. Failures are now contained per sample
  and summarised at the end.

### Added

- `src/utils/shellquote.py` — `shq()`, the single helper used to quote paths in shell
  command strings.
- Pre-flight input validation with `[ok]` / `[invalid]` / `[skip]` reporting, plus
  `--verify-gzip` to read gzipped inputs to EOF during validation.
- `[ok]` / `[empty]` / `[fail]` per-sample summary and a `STATUS: SUCCESS |
  PARTIAL_SUCCESS | FAILURE` line for the enrichment chain.
- `--EGGNOG_DATA_DIR` and `--kojson` on the `Onestep` command, so a local eggNOG
  database and KO hierarchy can be used without being downloaded.
- `data-index.env` — a single place that records which Release hosts
  `data-index.tar.gz` (see *Changed*).
- `CHANGELOG.md` (this file).

### Changed

- **Exit-code convention** shared by the enrichment chain: `0` success, `1` failure,
  `2` PARTIAL_SUCCESS (some sub-tasks failed, the rest produced results), `3` NO_TERMS
  (ran fine but found nothing significant — an empty result, not an error). A partial
  failure now always ends in a non-zero exit code *and* an explicit status line.
- Identification no longer aborts a batch when a sample legitimately yields no
  candidate; such samples are skipped, reported, and recorded as empty results.
- **The data index is decoupled from the software version.** `data-index.env` is now the
  single source of truth read by both workflows and by `test.sh`; the README and
  `CONTRIBUTING.md` quote the same tag in their download commands. When the index has to
  be updated, publish it as a new Release and change that one line — no software version
  bump and no changes to the workflows are needed.
- Dependency and environment files are unchanged in this release (`mirdp3_environment.yml`,
  `conda.recipe/meta.yaml`).
