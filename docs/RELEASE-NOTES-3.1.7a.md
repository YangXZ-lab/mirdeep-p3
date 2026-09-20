> Historical checkpoint: superseded by docs/RELEASE-NOTES-3.1.7b.md (same code, defaults flipped).

# Release notes: miRDeep-P3 3.1.7a

Status: **local only**, not pushed to GitHub (as instructed).
Baseline: /nfs/users/zhaojw/software/miRDeep-p3 (3.1.6d, never pushed either).
This working copy: /nfs/users/zhaojw/software/miRDeep-p3-3.1.7a.
Public repo: /nfs/users/zhaojw/software/github/git/mirdeep-p3 (still 3.1.6c).

## 1. What 3.1.7a adds

One feature, one engine: a psRNATarget-compatible target scanner, reachable as
`analysis Target_finder --engine psrna` and as `analysis Onestep --target_engine psrna`.
The historical engine is untouched in behaviour and remains the default.

Why it was needed: the historical engine applies psRNATarget scoring to the HSPs that
ssearch36 happens to report, and ssearch36 is called with `-b 200` (at most 200 HSPs per
miRNA). Measured against archived psRNATarget output on Arctium_lappa it recovers 352 of
562 pairs (62.6%). The new engine searches every position and recovers 559 (99.5%).

## 2. File-level diff against 3.1.6d

| File | Change | Lines |
|---|---|---|
| `src/version.py` | version string 3.1.6d -> 3.1.7a | 1 |
| `src/MiRNATarget/psrna_scan.py` | **new**: standalone psRNATarget-compatible scanner | 369 |
| `src/commands/analysis/target_finder.py` | dual engine, deterministic output order, `--engine/--allow-gap/--psrna-detail/--max-hsp` | 258 (was 158) |
| `src/commands/analysis/onestep.py` | `--target_engine` forwarding (one argument, one command line) | +5 |
| `docs/PSRNA-ENGINE.md` | **new**: algorithm, rules, validation, limits | 111 |
| `docs/RELEASE-NOTES-3.1.7a.md` | **new**: this file | - |
| `CHANGELOG.md` | 3.1.7a section (Added / Changed / Validation) | +67 |
| `README.md` | Target_finder section rewritten: psrna example, four new options, output note, notes 1-5 | +60 |

Untouched: identification, annotation, all other analysis subcommands, every R script,
data/, examples/, tests/.

## 3. Behaviour changes that can be observed

1. `target_finder.tsv` row order is now deterministic (mature FASTA order, then target,
   then target start). Before, two runs with identical input could differ in row order
   because ssearch36 `-T` output order is not stable. Content is unchanged; a user who
   sorted the file downstream sees nothing different.
2. `Target_finder` checks tools per engine: `ssearch` needs ssearch36 + seqkit, `psrna`
   needs numpy only. A machine without ssearch36 can now run target prediction.
3. `--GUs` is rejected as meaningless / ignored under `psrna` (psRNATarget V2 has no
   G:U cap), and `--max-hsp` is rejected under `psrna`. Wrong-engine options fail
   fast instead of being silently ignored.
4. `Onestep` accepts `--target_engine`; its default keeps the old path exactly.

## 4. How it was validated

* Identity with the research prototype: the CLI on Arctium_lappa returns the same 562
  site set as the stand-alone prototype developed in
  /nfs/users/zhaojw/software/psRNATarget_test (zero set difference, only row order).
* Against archived psRNATarget web output, four datasets:

| dataset | rows (psRNATarget) | rows (psrna) | missing | extra | identical, alignment excluded |
|---|---:|---:|---:|---:|---:|
| Arctium_lappa | 563 | 562 | 3 | 2 | 99.47% |
| Helianthus_annuus | 1,724 | 1,729 | 3 | 8 | 99.54% |
| Lactuca_sativa | 2,531 | 2,588 | 7 | 64 | 99.64% |
| Bidens_alba | 12,296 | 12,776 | 151 | 631 | 98.45% |

* Regression of the old engine: 3.1.6d and 3.1.7a produce identical `target_finder.tsv`
  after sorting, on the same input and settings.

## 5. Later changes

### 5.1 Short-term (recommended next)

* `test.sh` does not cover target prediction. Add a small fixture (a dozen miRNAs, a few
  hundred CDS) and assert row counts for both engines, so a future refactor cannot break
  the six-column contract silently.
* `analysis` subcommands still hard-code `Rscript~/`python~/`fimo~/`ssearch36~/`seqkit` paths
  (issue I-03 in docs/MIRDEEP-P3-ISSUES.md). Adding `--engine psrna` removed one such
  dependency; the rest remain.

### 5.2 Design questions left open

* Whether `psrna` should become the default engine. Keeping `ssearch` as default is the
  conservative choice (no silent change for existing pipelines), but the accuracy gap
  (62.6% vs 99.5%) argues for flipping it in a future minor release with a note.
* `--allow-gap` is off by default because enabling it adds sites psRNATarget does not
  report. If indel handling matters more than exact agreement, flip the default there.

## 6. Checklist for syncing to GitHub

The public repo is at 3.1.6c, so a sync carries **two** releases: 3.1.6d (expression
matrix fixes, already documented in CHANGELOG.md) and 3.1.7a (this one). From
/nfs/users/zhaojw/software/miRDeep-p3-3.1.7a copy only:

    src/version.py
    src/expression_matrix.py
    src/commands/annotation.py
    src/commands/analysis/onestep.py
    src/commands/analysis/target_finder.py
    src/MiRNATarget/psrna_scan.py        (new)
    docs/PSRNA-ENGINE.md                 (new)
    docs/RELEASE-NOTES-3.1.7a.md         (new)
    CHANGELOG.md
    README.md

Do **not** copy: `data/index/` (index tarball, not in git), `tests/data/` (49 MB fixture),
`TEST-3.1.6d.md` (local test log), any `__pycache__`. Then: create tag `v3.1.7a` (source)
and `v3.1.7a-full` (Docker/offline), and rebuild the Docker image, since the image tag in
README points at the released version.

Before pushing, run `./test.sh` from the 3.1.7a tree; note that it downloads the data
index into /tmp, which this account forbids -- set `TMPDIR` to a directory under
/nfs/users/zhaojw/Genburten/temp instead.
