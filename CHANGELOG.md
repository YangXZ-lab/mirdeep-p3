# Changelog

All notable changes to MirDeep-P3 are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Release tags follow the scheme `v<version>` (source, used by the conda recipe) and
`v<version>-full` (Docker images + offline bundle, see `.github/workflows/docker-release.yml`).
Releases before 3.1.6c are listed on the
[GitHub Releases page](https://github.com/YangXZ-lab/mirdeep-p3/releases).

## [Unreleased] - 3.1.7b

3.1.7b makes that second engine the default, turns on indel tolerance, and adds a
fixture that exercises the engine in test.sh.

Adds a second alignment engine to miRNA target prediction. The historical engine
(ssearch36 + the MiRNATarget parser) applies psRNATarget scoring rules, but only to the
HSPs ssearch36 reports, and ssearch36 -b 200 caps how many HSPs one miRNA can return. In
practice it agrees with psRNATarget far less than assumed: on Arctium_lappa it recovers
352 of 562 archived miRNA-target pairs (62.6%), while the new engine recovers 559 (99.5%).

### Added

- **--engine psrna for analysis Target_finder** (`src/commands/analysis/target_finder.py`).
  A built-in psRNATarget-compatible scanner that needs neither ssearch36 nor seqkit. It
  reproduces the psRNATarget V2 penalty scheme (1.5x seed mismatch, 0.5 G:U, 2.0 gap plus
  1.0 when the gap is on the miRNA side, only positions 1-19 scored), the Translation /
  Cleavage call, and the per-pair Multiplicity. Since 3.1.7b this is the default engine
  (see Changed); --engine ssearch keeps the old path.

- **src/MiRNATarget/psrna_scan.py** (new, 369 lines). The engine as a standalone script, so
  it can be run outside miRDeep-P3: -m mature.fa -t targets.fa -o out.tsv -e 2.5 -p 32.
  It writes the native 14-column psRNATarget table. Only numpy is required.

- **--allow-gap / --no-allow-gap** (psrna only). Whether to report sites carrying a
  single-base indel. On by default since 3.1.7b: psRNATarget reports such sites, so
  excluding them cost recall, and keeping them costs only a few percent of extra rows
  (measured below). --no-allow-gap restores the indel-free subset.

- **--psrna-detail** (psrna only). Additionally write the native 14-column table as
  <output>/target_finder.psrna.tsv, for users who want the psRNATarget fields (UPE,
  inhibition, multiplicity, aligned fragments) that the six-column file drops.

- **--max-hsp** (ssearch only). The ssearch36 -b value, default 200 (historical
  behaviour). Raising it reduces the truncation described above.

- **--target_engine {ssearch,psrna}** on analysis Onestep, forwarded to Target_finder.

- **docs/PSRNA-ENGINE.md**: the rules as implemented, the search strategy, the four-dataset
  validation table, and where the residual differences come from.

- **Target-prediction fixture in the test suite** (tests/data/target_mini_mature.fa,
  tests/data/target_mini_cds.fa, section 9 of test.sh). 20 mature miRNAs and 300 CDS
  extracted from the Arctium_lappa validation set, small enough for CI. test.sh asserts
  the exact row count (120), the six-column shape, byte-identical repeat runs, that
  --no-allow-gap leaves this fixture unchanged, and (when ssearch36 is installed) that
  the historical engine still returns at least 100 rows.

### Changed

- **analysis Target_finder now defaults to --engine psrna.** This is a behaviour
  change: the default no longer needs ssearch36 or seqkit, and the recall is much higher
  (see Validation). Pass --engine ssearch to keep the 3.1.6x behaviour. Onestep got the
  same default through --target_engine.

- **--allow-gap is on by default** (psrna engine only). Sites with a single-base indel
  are scored the way psRNATarget scores them; --no-allow-gap excludes them. On the
  validation datasets this adds 0.0-4.3% rows and recovers missed pairs (see Validation).

- **target_finder.tsv is now written in a deterministic order**: mature FASTA order, then
  target ID, then target start. Previously the row order changed between runs with the same
  input because ssearch36 -T output order is not stable (two runs of 3.1.6d on the same
  data differed in row order but not content). Downstream consumers that sorted the file
  themselves see no change; consumers that relied on the incidental order should re-check.

- **Target_finder only requires seqkit for --engine ssearch.** The psrna engine needs
  numpy only, and the tool check now reflects the selected engine.

- **--GUs is documented and enforced as ssearch only.** psRNATarget V2 has no G:U cap;
  passing a non-default value with --engine psrna prints a note and ignores it.

### Validation

Four plant datasets, comparing against archived psRNATarget web output. "Identical"
compares every field except the alignment string, which psRNATarget renders only over the
HSP extent.

With the 3.1.7b defaults (psrna engine, --allow-gap on):

| dataset | psRNATarget rows | psrna rows | delta | missing | extra | identical (excluding alignment) |
|---|---:|---:|---:|---:|---:|---:|
| Arctium_lappa | 563 | 579 | +16 (+2.8%) | 1 | 17 | 99.82% |
| Helianthus_annuus | 1,724 | 1,773 | +49 (+2.8%) | 1 | 50 | 99.65% |
| Lactuca_sativa | 2,531 | 2,639 | +108 (+4.3%) | 5 | 113 | 99.72% |
| Bidens_alba (first 120 miRNAs) | 3,269 | 3,270 | +1 | 42 | 43 | 98.35% |

Bidens_alba was compared on a 120-miRNA subset over all 165,431 CDS, not on its full
511 miRNAs, so the indel-free run below is quoted for that dataset for reference.

With --allow-gap off (the 3.1.7a behaviour):

| dataset | psRNATarget rows | psrna rows | missing | extra | identical (excluding alignment) |
|---|---:|---:|---:|---:|---:|
| Arctium_lappa | 563 | 562 | 3 | 2 | 99.47% |
| Helianthus_annuus | 1,724 | 1,729 | 3 | 8 | 99.54% |
| Lactuca_sativa | 2,531 | 2,588 | 7 | 64 | 99.64% |
| Bidens_alba | 12,296 | 12,776 | 151 | 631 | 98.45% |

The ssearch engine was regression-tested against 3.1.6d: with the new ordering applied,
the two outputs are byte-identical, so the psrna work did not change the historical path.

## [Unreleased] - 3.1.6d

This release fixes two silent losses in the expression matrices: miRNAs that were
annotated but never quantified, and reads that were quantified only partially.

### Fixed

- **Annotated miRNAs silently missing from the expression matrices**
  (`src/expression_matrix.py`). The row set of `*-mature.count` / `*-mature.exp` was
  derived from the alignments alone, so a miRNA without a single qualifying read got no
  row at all, while remaining in `*-basic-info`. The two files could therefore not be
  joined. On a 21-sample rice dataset this produced 446 annotated miRNAs against 427
  matrix rows, with no warning anywhere.

  `expression_matrix.py` now accepts `-fa <prefix>-basic-info-mature-1nt.fa`. When given,
  every reference miRNA is written to the matrix (zero-filled when no read passes the
  filter) and the row order follows the reference FASTA, which is the same order as
  `*-basic-info`. `mirdeep-p3 annotation` passes this argument automatically. miRNAs that
  already had reads keep their previous values unchanged.

- **Over-strict coverage filter** (`src/expression_matrix.py`). Any read ending more than
  2 bp short of the mature reference 3' end was discarded, which removed genuine
  3'-truncated reads; plant sRNA-seq libraries are dominated by 19-21 nt reads, so one
  missing base was enough to drop a miRNA entirely. Both ends are now configurable through
  `--max-5p-trunc` (default 2, unchanged) and `--max-3p-trunc` (default 3, previously
  hard-coded to 2), exposed by `expression_matrix.py` and `mirdeep-p3 annotation`.
  `--max-3p-trunc 2` reproduces the <=3.1.6c behaviour exactly.

  Measured on the same dataset (21 samples, 22,325 alignment lines):

  | setting | alignments kept | rows | non-zero rows | total counts |
  |---|---|---:|---:|---:|
  | 5'<=2, 3'<=2 (3.1.6c) | 66.6% | 427 | 427 | 46,444,745 |
  | 5'<=2, 3'<=2, with `-fa` | 66.6% | 446 | 427 | 46,444,745 |
  | 5'<=2, 3'<=3 (new default) | 78.7% | 446 | 431 | 46,831,031 |
  | 5'<=2, 3'<=4 | 83.2% | 446 | 432 | 49,781,809 |
  | 5'<=4, 3'<=4 | 95.8% | 446 | 433 | 50,359,058 |

  3'<=4 raises the MIR167 family by 41% and loosening the 5' end to 4 bp raises
  MIR12577a-h by 61x, because 5'-truncated reads then get recruited from unrelated loci.
  Neither is recommended, hence the default of 3.

- **`--min-expr` not reachable from `Onestep`** (`src/commands/analysis/onestep.py`).
  `Differential_expression` filters the matrices to miRNAs reaching `--min-expr` RPM in at
  least one sample (default 5.0) before PCA, correlation and DESeq2, but `Onestep` called
  that subcommand without the flag, so the threshold was stuck at 5.0 for the whole
  pipeline. `Onestep` now accepts `--min-expr` and forwards it; a negative value is
  rejected before the output directory is cleared and before any step runs.

### Changed

- `expression_matrix.py` reports the number of alignments kept and dropped, plus the
  number of zero-filled rows, on stderr. The coverage filter is no longer silent.

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
