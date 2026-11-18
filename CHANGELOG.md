# Changelog

All notable changes to MirDeep-P3 are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Release tags follow the scheme `v<version>` (source, used by the conda recipe) and
`v<version>-full` (Docker images + offline bundle, see `.github/workflows/docker-release.yml`).
Releases before 3.1.6c are listed on the
[GitHub Releases page](https://github.com/YangXZ-lab/mirdeep-p3/releases).

## [Unreleased] - 3.1.8d

3.1.8d closes the loop that 3.1.8c left open: a corrected mature can end up one
base outside its own precursor, and until now that was reported but not fixed.

### Fixed

- **A corrected mature can overhang the annotated precursor**
  (`src/utils/precursor_repair.py`, `src/repair_precursor.py`,
  `src/commands/annotation.py`). The precursor is excised flush with the mature
  it was built around; move the mature one base downstream and its 3' end lands
  one base outside the hairpin. Nothing crashes, but the mature stops being a
  substring of the precursor, so the stored dot-bracket describes a window that
  no longer contains the mature it annotates, and the star is placed relative to
  a frame that is now off by one.

  Measured on the rice dataset: 0 of 446 loci overhang before the correction,
  2 after (`Nip-MIR166a`, `Nip-MIR166b`), both by exactly 1 nt. Of 1420 loci
  across four species, only those two are affected.

  The new step runs only when `--fix-mature-boundary` actually rewrote the
  table. For each overhanging locus it widens the precursor on the short side
  (bases read from the genome), re-folds the widened precursor and its
  +/-20 nt window with RNAfold, and recomputes the star. The widening is
  verified to have only *added* bases -- the old precursor must still be present
  at the offset the coordinates imply, otherwise the genome is a different
  assembly and the row is left alone.

  The mature sequence is untouched, so the bowtie reference and the alignments
  stay valid; no third round of alignment is needed. On the rice dataset the
  repair changes 2 of 446 rows and only columns 6-13 (the two windows); the
  other 444 rows stay byte-identical and the star does not move, because the
  extra base goes on the transcript 3' arm and the star sits on the other one.

### Added

- **`src/utils/genome_fasta.py`** -- random access to a genome FASTA with a
  self-built byte-offset index, cached beside the file (falling back to a temp
  directory when that is read-only). `samtools faidx` and `bedtools getfasta`
  both want to write an index next to the FASTA, which is not always allowed;
  pyfaidx would be a new dependency. A 360 MB genome costs one scan (~2 s) and
  a few hundred bytes per query.
- **`src/utils/precursor_repair.py`** -- the widening logic, stdlib only:
  the interval arithmetic for "which side is short", a transcript-offset check
  that the old precursor still sits where it should, the guard rail
  (`--max-extension`, default 5 nt per side), and star recomputation.
- **`src/repair_precursor.py`** -- CLI, also usable on its own against any
  22-column basic-info plus its genome. Writes `<prefix>-precursor-repair.tsv`
  with the old/new bounds, how many bases were added on each side, whether the
  star moved, and how many columns changed.
- **`tests/test_precursor_repair.py`** -- stdlib-only unit tests: the genome
  fetcher against a synthetic multi-line FASTA on both strands (an off-by-one in
  the byte arithmetic is invisible on a single-line record), the widening
  decision for the high/low/both-sides cases, the refusal paths (too long,
  sequence mismatch, RNAfold failure), and a check that writing a repair touches
  only the window and star columns.

### Changed

- **`annotation --fix-mature-boundary`** now also repairs precursors, keeping
  the pre-repair table as `<prefix>-basic-info.before-precursor-repair` and
  writing a second report, `<prefix>-precursor-repair.tsv`.
- **The basic-info interval is documented as 0-based.** 3.1.8c called it
  1-based. Nothing computed wrong -- only differences of two coordinates are
  used -- but the genome check settles it: the stored precursor and mature
  sequences reproduce IRGSP-1.0 exactly for 446/446 loci only under the 0-based
  reading (19/19 on the Arabidopsis test data). Anyone slicing a genome with
  these coordinates, as the new repair step does, has to use 0-based.

### Notes

- Default behaviour is unchanged: with no `--fix-mature-boundary` the table is
  byte-identical to 3.1.8b/3.1.8c (`f9223f3c...` on the mini dataset) and no
  repair report is produced.
- The repair never touches the mature interval, the names, the family, the
  MIR_ID or the cluster-relevant columns.
- Four species re-checked for scope: a correction changes at most 2 rows of 446,
  and only the mature triplet; the cluster file and the other 444 quantified
  loci are unchanged.

## [Unreleased] - 3.1.8c

This release fixes a silent 1-nt error in the mature 5' boundary.

### Fixed

- **Mature 5' boundary can be one base off** (`src/utils/mature_boundary.py`,
  `src/audit_mature_boundary.py`). The mature interval comes from the prediction
  and is then adjusted per arm by `utils/coords.py`; that adjustment can leave it
  one base longer or shorter at the 5' end than the reads that support the locus.
  Nothing fails when this happens, and the expression matrix absorbs it as a
  1-nt offset rather than reporting it.

  Measured on four annotated datasets (Arabidopsis 135, Gnetum 513, rice 446,
  Scaevola 344 records; 1420 loci with usable coverage), 47 loci (3.3%) have a
  dominant read offset that is not canonical. Three of them are confirmed
  annotation errors against miRBase: `Nip-MIR166a`, `Nip-MIR166b` (annotation
  keeps one extra base at the 5' end, the reads are the miRBase mature) and
  `Gmo-MIR159a`. `Nip-MIR166a/b` are the two most abundant members of that
  family, 1.6M counts each, so the error is not numerically small.

  What annoys the diagnosis: 24 of the 47 loci are real 5'-truncated molecules,
  not annotation errors — their annotated mature already sits on the arm
  boundary, so there is no upstream room for the annotation to have drifted.
  The audit therefore reports rather than assumes.

### Added

- **`src/utils/mature_boundary.py`** — the coordinate arithmetic for this
  problem, stdlib only: precursor-relative <-> genomic conversion honouring
  strand, dominant-offset selection with a deterministic tie-break towards
  canonical, and a guarded correction that refuses anything leaving the
  precursor or producing a mature outside 18-26 nt.
- **`src/audit_mature_boundary.py`** — CLI. Reads `*.mature.aln`, keeps only
  reads that reach the 3' end (`--max-3p-gap`, default 2) because a read
  truncated at its 3' end says nothing about where the 5' starts, and writes
  `<prefix>-boundary-audit.tsv` with the dominant offset, its share, an
  arm-placement label and the proposed interval for every locus.
- **`tests/test_mature_boundary.py`** — stdlib-only unit tests: relative/genomic
  round trip on both strands, the trim convention, tie-breaking, the correction
  guards, and the strand mirroring that a wrong implementation gets backwards
  without failing.
- **`annotation --fix-mature-boundary`** — opt-in correction. When set, the
  audit runs with `--apply`, the pre-correction table is kept as
  `<prefix>-basic-info.before-boundary-fix`, the reference and alignments are
  rebuilt from the corrected table, and a second audit
  (`<prefix>-boundary-audit.after-fix.tsv`) verifies the offsets are now
  canonical.
- **`annotation --boundary-min-fraction` / `--boundary-min-weight`** — the
  evidence thresholds, default 0.60 and 10.0 read weight.

### Changed

- **The annotation part 2 order** (`src/commands/annotation.py`): align ->
  audit (+ optional correct) -> cluster -> expression. The cluster used to run
  before quantification; it now runs after the audit so that a corrected table
  reaches every downstream product instead of only the matrices.

### Notes

- The audit runs on every annotation, correction does not. Defaults leave the
  table byte-identical to 3.1.8b; only `--fix-mature-boundary` rewrites it.
- Only +/-1 nt corrections are applied. A locus that is off by two or more is
  reported as `shift_beyond_limit` and left alone, because at that distance the
  dominant offset is more likely to be 5'-degradation than a boundary slip.
- The audit measures two windows, and they are not interchangeable. Candidates
  are built in the precursor +/-20 nt window, which is wide enough that a mature
  shifted by one base still lands inside it; the arm-placement label and the
  containment check use the real hairpin (mature + loop + star). A 22 column
  table carries both. Reading only the wide window made the arm label a
  20 nt-too-large yardstick and hid the case below.
- A correction can push the mature one base past the annotated precursor end,
  which is what happens when the hairpin was excised flush with the old (wrong)
  3' end: `Nip-MIR166a` and `Nip-MIR166b` on rice are exactly that. They stay
  corrected, because the reads are the evidence, and the run log names them
  under `corrected mature overhangs the annotated precursor: N (...)`. Extending
  the precursor would also mean re-folding the hairpin and rewriting the stored
  structure, so this is reported as a known limitation instead of edited
  silently.
- A correction touches the mature interval and nothing else: measured on the
  rice dataset, 2 of 446 rows change and only columns 16/17/18 move, the other
  444 rows stay byte-identical, the cluster file is unchanged, the expression
  reference differs in exactly those 2 entries, and re-quantifying the other 444
  loci returns the same counts.

## [Unreleased] - 3.1.8b

3.1.8b is housekeeping: the shared coordinate arithmetic gets a home and a test,
the combined `--step` help is generated from the real parsers, the tree gains a
standard console entry point, and the installed data index finally says which
release it is. No algorithm, no default and no output format changed.

### Added

- **`src/utils/coords.py`** (new, O-05). `compute_star()` and
  `map_prediction()` now hold the coordinate arithmetic that used to be inlined
  in `src/primary-basic-info.py` and `src/mod_rm_redundant_meet_plant.py`
  (116 lines removed from the two scripts). `map_prediction()` raises
  `ValueError("negative or zero positions")`; the caller keeps its original
  warning text and skips that record exactly as before.
- **`tests/test_coordinates.py`** (new, standard library only). Covers
  `compute_star()` on both strands and both mature placements, the regression
  that a star interval must be as long as its sequence (the 3.1.6c bug), and
  `map_prediction()` for the four strand x arm combinations against
  hand-computed coordinates.
- **`src/utils/data_index.py`** (new, O-10). Reports the installed data-index
  release and compares it with `data-index.env`; the states are `missing`,
  `unknown`, `unrecorded`, `match` and `mismatch`. `identification` and
  `annotation` print one `[index] data index: ...` line at start-up, the batch
  manifest records the whole state, and a changed tag invalidates `--resume`.
  `tests/test_data_index.py` pins every state.
- **`pyproject.toml` + `mirdeep_p3/` package** (new, O-09, conservative). A
  standard console script (`mirdeep-p3 = "mirdeep_p3.cli:main"`) that locates
  the checkout and runs the existing entry point, so `pip install -e .` works.
  The tree layout, `data/` and `scripts/` are untouched and nothing is packaged
  into a wheel.
- **`MIRDEEP_P3_HOME`** environment variable, honoured by `./mirdeep-p3` too,
  for installations that keep the tree outside the package directory.

### Changed

- **Combined `--step` help is generated** (O-06) from the selected steps own
  parsers. `--step 1,3`, `--step 2,3` and `--step 1,2,3` now print a real option
  list instead of the "(Full combined help will be implemented later.)"
  placeholder. The 55-line hand-written 1+2 help is gone; the generated output
  groups options by the steps that accept them and marks options whose help text
  differs between steps.
- `src/commands/analysis/__init__.py` exposes `build_parser()` so the combined
  help reuses the real parser; `analysis.main()` behaviour is unchanged.
- `test.sh` gained a CLI/help/unit section (exit code 7) and records
  `data/index/.index-version` when it unpacks the index.

### Fixed

- The first generated-help build crashed on `--step 1,3`, `--step 2,3` and
  `--step 1,2,3` with `AttributeError: ArgumentParser object has no attribute
  help`; the subcommand help now comes from `_choices_actions` and the aliases
  (`DE`, `FA`) are collapsed.

### Documentation

- **`docs/DATA-INDEX.md`** (new). The data index is versioned on its own axis, so
  the document spells out the two axes, the `.index-version` marker, the five
  states, every file that reads `DATA_INDEX_RELEASE`, and an 8-step checklist for
  publishing a new index (release asset, `data-index.env`, the two `wget`
  commands, the marker, `test.sh`, `CHANGELOG.md`, the user-facing note).
- **`README.md`** now opens with a *Current versions* block quoting the software
  version and the current data-index tag, linking `data-index.env` (single source
  of truth) and `docs/DATA-INDEX.md`.
- **`test.sh` section 11** guards those docs: it fails with exit code 7 when the
  README block disagrees with `src/version.py` or with `data-index.env`, or when
  `docs/DATA-INDEX.md` no longer mentions the current tag. Verified by running the
  guard against a copy of the tree with a wrong version and with a wrong tag.

### Validation

| Check | Result |
|---|---|
| `python3 tests/test_coordinates.py` | all checks passed |
| `python3 tests/test_data_index.py` | 5/5 passed |
| `--step 1,2` / `1,3` / `2,3` / `1,2,3` `-h` | all exit 0; `1,3` has no `--prefix_miRNA`, `2,3` does, `1,2,3` prints all three groups |
| identification on the mini dataset (12 stages) | every output md5 identical to 3.1.8a (`mini_nr_predictions` 3ed474325e0c0a3ff58f0ef738106229, `mini_filter_P_prediction` 2e63d35040c892d3fc67797ba2bab877) |
| annotation on the same data | `mirdp3-test-basic-info` f9223f3c2b120d62ec24fa60095aca09, `-basic-info-cluster` db799e4e78ea901bcb75fb2d0f160dbd — identical to 3.1.8a |
| index reporting | `[index] data index: ... (no .index-version; data-index.env expects mirdeep-p3-v3.1.4c-full)`; batch manifest carries `params.data_index`, sample manifest carries the installed tag |
| `./test.sh -t 12 -k` | ALL TESTS PASSED (identification, annotation, Target_finder psrna + ssearch, Stat, DE, CLI/help/unit) |

## [Unreleased] - 3.1.8a

3.1.8a turns the identification pipeline into a declared list of stages and gives
every run a machine-readable record. The science is unchanged: on the mini dataset
the five key outputs are byte-identical to 3.1.7c (see Validation).

### Added

- **Declarative stage runner** (`src/commands/stage_runner.py`, new). The fifteen
  inline steps of `identification` are now twelve `Stage` objects (id, name,
  artifacts, run). Each stage is timed, its artifacts are recorded in the manifest,
  and "empty result" is expressed as `raise StageEmpty(reason)` instead of a special
  return value, so the `ok / empty / failed` semantics of 3.1.7c are preserved
  verbatim (the six original reason strings are unchanged).

- **`--resume` for identification.** A stage is skipped when the recorded manifest
  still matches (same version, parameters, genome, input fingerprint and index) and
  every artifact of that stage is present and non-empty. Changed inputs invalidate
  the manifest and the sample runs from the start. `--clean` is ignored while
  resuming, because `temp/` holds the intermediates resuming depends on.

- **Run manifest.** `<out>/<prefix>/<prefix>.run.json` records schema version,
  software version, parameters, tool versions (`bowtie`, `bowtie-build`, `RNAfold`,
  `samtools`, `trim_galore`), input fingerprint (size, mtime, md5), genome
  fingerprint, and per-stage status, artifacts and duration.
  `<out>/mirdp3-identification-<time>.manifest.json` summarises the batch (one entry
  per sample). Both are written atomically (`.part` + `os.replace`).

- **`analysis` smoke test in `test.sh`**: `Stat` on the annotation output plus
  `Differential_expression` on a bundled 180-miRNA, six-sample fixture
  (`tests/data/analysis_mini.count`, `tests/data/analysis_mini.exp`). The DE part is
  invoked through its `DE` alias on purpose, so the alias fix below stays covered.
  New exit code `6`.

### Fixed

- **`analysis DE` / `analysis FA` aliases did nothing**
  (`src/commands/analysis/__init__.py`). The subparsers are registered with aliases,
  but dispatch compared against the full name, so invoking an alias fell through to
  the help branch and exited 1 (no error message, no results). The command name is
  now normalised before dispatch.

### Changed

- `identification --clean` is silently disabled when `--resume` is used, with an
  explicit log line, instead of deleting the tree that makes resuming possible.

- The batch manifest is named with second resolution
  (`mirdp3-identification-<mmddyy-HHMMSS>.manifest.json`) so two runs started in the
  same minute no longer overwrite each other's record. The `.pipe` name keeps its
  minute-resolution form for compatibility with existing downstream scripts.

### Validation

Byte-identical outputs on the mini dataset (`tests/data/mini.fq` + `mini_genome.fasta`,
`-t 12 -p 4`), comparing 3.1.7c and 3.1.8a:

| Output | md5 |
|---|---|
| `mini_nr_predictions` | `3ed474325e0c0a3ff58f0ef738106229` |
| `mini_nr_predictions.bed` | `e17fc655d0b05419a1ec21b16a14b996` |
| `mini_filter_P_prediction` | `2e63d35040c892d3fc67797ba2bab877` |
| `mini.processed.fa` | `760c9775cc6acdf82bdc4ef4184743c7` |
| `mini.total_reads` | `2f3610f0470b89e45076bbfcc89dd538` |

Resume behaviour (mini dataset):

| Scenario | Result |
|---|---|
| `--resume` on an unchanged run | 12/12 stages `[resume] ... skip`, seconds, outputs unchanged |
| `--resume` after deleting `mini_nr_predictions.bed` | 11 skipped, `4.15` re-ran, rebuilt bed md5 matches |
| edited manifest signature | `input or parameters changed, not resuming`, full re-run |

`./test.sh -t 12 -k` (with `Rscript` present): identification, annotation, target
prediction (psrna 120 rows deterministic, ssearch 112 rows), `Stat` (4 figures) and
`Differential_expression` (180 miRNAs, DEG lists written) — `ALL TESTS PASSED`, exit 0.

## [Unreleased] - 3.1.7c

3.1.7c is a maintenance release: it closes the entries of the issue ledger that were
approved for this cycle, makes the analysis commands fail early with an actionable
message when an external tool is unusable, and stops leaving large temp directories
behind.

### Added

- **External tool paths for every analysis subcommand** (`src/utils/toolset.py`, new).
  Stat, TFBS, Target_finder, Differential_expression, Functional_analysis and Onestep
  accept `--rscript`, `--fimo`, `--bedtools`, `--ssearch36`, `--seqkit`,
  `--emapper` and `--python`. Each command prints an `External dependency check:`
  block with the resolved paths and exits before doing any work if a required tool is
  missing or not executable. Onestep forwards whatever it was given to its sub-steps.

- **`--keep-temp` for Target_finder** (`src/commands/analysis/target_finder.py`).
  The temp directory is removed when the run finishes; `--keep-temp` keeps the
  intermediate alignments, matching the behaviour of the other analysis commands.

- **`--clean-index` for identification** (`src/commands/identification.py`).
  `--clean` now really removes only `<prefix>/temp`, and the bowtie index is deleted by
  the new flag. When the index was supplied with `-d` it is left in place and the skip
  is logged.

- **`src/miRNA_consistency.py`** (restored). The script behind
  `annotation --consistency` and `--common` was missing from the repository, so both
  flags aborted with `FileNotFoundError`.

- **`scripts/dedup_score_pair.py`** (restored). Referenced by
  `docs/core-dataset.md` but absent, which broke the documented core-dataset
  evaluation pipeline.

### Changed

- **`src/check_python_deps.py` now reports what it actually checks.** The pipeline
  core needs no third-party package; numpy is reported as an optional dependency of
  `Target_finder --engine psrna`. The previous list contained only standard-library
  modules, so the check could never fail.

- **`--step` stripping in `./mirdeep-p3`.** The two `del cli_args[...]` calls assumed
  a fixed offset; the flag and its value are now filtered by index, so `-h --step 1,2`
  and `--step 1,2 -h` produce identical help.

### Fixed

- Executable bits: `mirdeep-p3`, `src/miRNA_consistency.py` and
  `scripts/dedup_score_pair.py` are committed as `100755`.

### Validation

- `./test.sh -t 12 -k` — ALL TESTS PASSED (identification, annotation, both target
  engines, determinism).
- Functional smoke suite (33 checks): help of the six analysis subcommands,
  `--step` parsing equivalence, both `miRNA_consistency` modes on a real basic-info
  file, Stat, TFBS on a BED input, Target_finder with and without `--keep-temp`
  (identical output), the ssearch engine, early abort on a bad `--rscript` /
  `--emapper`, and DE on a negative-binomial 3v3 matrix (237 miRNAs tested,
  33 up / 22 down, temp cleaned).


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
