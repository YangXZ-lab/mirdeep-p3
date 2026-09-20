# Release notes: miRDeep-P3 3.1.7b

Status: **local only**, not pushed to GitHub (as instructed).
Baseline: /nfs/users/zhaojw/software/miRDeep-p3 (3.1.6d, never pushed either).
This working copy: /nfs/users/zhaojw/software/miRDeep-p3-3.1.7a (same tree, version bumped in place as instructed).
Public repo: /nfs/users/zhaojw/software/github/git/mirdeep-p3 (still 3.1.6c).
Supersedes the checkpoint note docs/RELEASE-NOTES-3.1.7a.md.

## 1. What 3.1.7b changes

One theme: the psRNATarget-compatible engine is now what users get by default.

1. `analysis Target_finder` defaults to `--engine psrna` (3.1.7a default: `ssearch`).
   `analysis Onestep` defaults to `--target_engine psrna` for the same reason.
2. `--allow-gap` is **on** by default (3.1.7a: off). `--no-allow-gap` restores the
   indel-free subset. The flag is now a `BooleanOptionalAction`, so both spellings
   exist and the default is visible in `-h`.
3. `test.sh` covers target prediction. Fixture: `tests/data/target_mini_mature.fa`
   (20 mature miRNAs) and `tests/data/target_mini_cds.fa` (300 CDS), both extracted from
   the Arctium_lappa validation set. The test asserts the exact row count (120), the
   six-column shape, byte-identical repeat runs, that `--no-allow-gap` leaves this
   fixture unchanged, and, when ssearch36 is installed, that the historical engine still
   returns at least 100 rows. New exit code 5 = target prediction failed.

## 2. Why the default engine flipped

Against archived psRNATarget web output, the historical engine recovers a minority of the
reported pairs because it can only score the HSPs ssearch36 returns (`-b 200` caps them):

| dataset | psRNATarget pairs | ssearch engine | psrna engine (indel-free) |
|---|---:|---:|---:|
| Arctium_lappa | 563 | 605 rows, 352 shared (62.6%) | 562 rows, 559 shared (99.5%) |

Keeping the weaker engine as the default only protects pipelines that would rather have
fewer, truncation-limited predictions. That is what `--engine ssearch` is for.

## 3. Why --allow-gap flipped

psRNATarget itself reports sites carrying a single-base indel. Excluding them (3.1.7a
default) was the larger source of missing pairs; including them costs a bounded number of
extra rows, because the search then also finds indel sites psRNATarget's HSP-based search
never reports.

| dataset | psRNATarget rows | gap off | gap on | missing off -> on | extra off -> on | identical off -> on |
|---|---:|---:|---:|---|---|---|
| Arctium_lappa | 563 | 562 | 579 (+2.8%) | 3 -> 1 | 2 -> 17 | 99.47% -> 99.82% |
| Helianthus_annuus | 1,724 | 1,729 | 1,773 (+2.8%) | 3 -> 1 | 8 -> 50 | 99.54% -> 99.65% |
| Lactuca_sativa | 2,531 | 2,588 | 2,639 (+4.3%) | 7 -> 5 | 64 -> 113 | 99.64% -> 99.72% |
| Bidens_alba (120 miRNAs, 165,431 CDS) | 3,269 | 3,263 | 3,270 (+0.2%) | 45 -> 42 | 39 -> 43 | 98.35% -> 98.35% |

Cost is bounded at a few percent of extra rows; recall improves or stays equal on every
dataset measured. Bidens_alba was compared on a 120-miRNA subset over all CDS, not on its
full 511 miRNAs.

## 4. File-level diff against 3.1.6d

| File | Change | Lines |
|---|---|---|
| `src/version.py` | 3.1.6d -> 3.1.7b | 1 |
| `src/MiRNATarget/psrna_scan.py` | **new**: standalone psRNATarget-compatible scanner (369 lines, numpy only) | 369 |
| `src/commands/analysis/target_finder.py` | dual engine, `DEFAULT_ENGINE = psrna`, `--engine/--allow-gap/--no-allow-gap/--psrna-detail/--max-hsp`, deterministic order | 299 (was 158) |
| `src/commands/analysis/onestep.py` | `--target_engine` forwarding, default `psrna` | +7 |
| `tests/data/target_mini_mature.fa` | **new**: 20 mature miRNAs (688 B) | 40 |
| `tests/data/target_mini_cds.fa` | **new**: 300 CDS (403 KB) | 600 |
| `test.sh` | section 9: target prediction, default engine + `--no-allow-gap` + ssearch | +55 |
| `docs/PSRNA-ENGINE.md` | **new**: algorithm, rules, validation, limits | 113 |
| `docs/RELEASE-NOTES-3.1.7b.md` | **new**: this file | - |
| `CHANGELOG.md` | 3.1.7b section (Added / Changed / Validation) | +80 |
| `README.md` | Target_finder section: psrna as the default example, `--no-allow-gap`, notes rewritten | +25 |

Untouched: identification, annotation, all other analysis subcommands, every R script,
`data/`, `examples/`.

## 5. Behaviour changes that can be observed

1. `target_finder.tsv` gains indel sites by default and loses them with `--no-allow-gap`.
   Row sets therefore differ from 3.1.6d even with identical inputs.
2. No ssearch36 / seqkit needed for the default path; only numpy.
3. Wrong-engine options fail fast: `--max-hsp` under psrna, `--allow-gap`/`--no-allow-gap`
   /`--psrna-detail` under ssearch.
4. `target_finder.tsv` row order is deterministic (mature FASTA order, then target, then
   start) for both engines.
5. `--GUs` is still accepted for `ssearch` and ignored with a note under `psrna`.

## 6. How it was validated

* Identity with the research prototype: on Arctium_lappa the psrna engine returns the
  same site set as /nfs/users/zhaojw/software/psRNATarget_test (zero set difference).
* Against archived psRNATarget output: the four tables above; scoring rules verified by
  replaying the archived alignment strings (563/563 and 1723/1724 Expectations
  reproduced) and by a 1920-combination grid search that picked this parameter set as the
  unique optimum.
* Determinism: two identical psrna runs are byte-identical (asserted in `test.sh`).
* Fixture: default engine 120 rows; `--no-allow-gap` 120; `--engine ssearch` 112.
* `ssearch` engine regression: 3.1.6d and the new tree produce the same rows after
  sorting, so the historical path is unchanged.
* `bash -n test.sh` passes. The full `test.sh` was not run end to end in this tree
  (identification needs `data/index`, and the script downloads it into /tmp, which this
  account forbids; set TMPDIR to a directory under /nfs/users/zhaojw/Genburten/temp).

## 7. Open items

* `analysis` subcommands still hard-code `Rscript` / `python` / `fimo` / `ssearch36` /
  `seqkit` paths (issue I-03 in docs/MIRDEEP-P3-ISSUES.md).
* `data/index/` in this tree is a local copy, not versioned; a GitHub sync must fetch the
  data-index Release instead.
* Docker image tag and README install instructions still point at 3.1.6c.

## 8. Checklist for syncing to GitHub

The public repo is at 3.1.6c, so a sync carries **two** releases: 3.1.6d (expression
matrix fixes, already in CHANGELOG.md) and 3.1.7b. From
/nfs/users/zhaojw/software/miRDeep-p3-3.1.7a copy only:

    src/version.py
    src/expression_matrix.py
    src/commands/annotation.py
    src/commands/analysis/onestep.py
    src/commands/analysis/target_finder.py
    src/MiRNATarget/psrna_scan.py        (new)
    tests/data/target_mini_mature.fa     (new, 688 B)
    tests/data/target_mini_cds.fa        (new, 403 KB)
    test.sh
    docs/PSRNA-ENGINE.md                 (new)
    docs/RELEASE-NOTES-3.1.7b.md         (new)
    CHANGELOG.md
    README.md

Do **not** copy: `data/index/` (index tarball, not in git), `tests/data/mini.fq` and
`tests/data/mini_genome.fasta` (49 MB / 31 MB fixtures already in the repo),
`TEST-3.1.6d.md` (local test log), any `__pycache__`.

Then: create tag `v3.1.7b` (source) and `v3.1.7b-full` (Docker/offline), rebuild the
Docker image, and update the image tag in README. Run `./test.sh` from the synced tree
first (with TMPDIR set as described above).
