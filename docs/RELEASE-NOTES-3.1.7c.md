# Release notes — miRDeep-P3 3.1.7c

Release date: 2026-11-18
Base: 3.1.7b (commit `1876aa1`)
Scope: maintenance — closes eight entries of the issue ledger (`docs/MIRDEEP-P3-ISSUES.md`).

## 1. What changed

| Ledger | Change |
|---|---|
| I-01 | `src/miRNA_consistency.py` is shipped again, so `annotation --consistency` / `--common` work |
| I-02 | `Target_finder` removes its temp directory by default; `--keep-temp` keeps it |
| I-03 | Every analysis subcommand resolves external tools through one helper and accepts explicit paths |
| I-04 | `check_python_deps.py` reports what it really checks (numpy is optional) |
| I-05 | `--step` is stripped from argv by index, not by a fixed offset |
| I-07 | `--clean` only removes `temp`; the new `--clean-index` removes the bowtie index |
| I-10 | `scripts/dedup_score_pair.py` is part of the repository again |
| I-13 | The entry point and the two restored scripts are committed as `100755` |

## 2. External tool paths

`analysis` subcommands no longer hard-code `Rscript`, `fimo`, `bedtools`,
`ssearch36`, `seqkit` or `emapper.py`:

```bash
mirdeep-p3 analysis Functional_analysis \
  --emapper /opt/eggnog/emapper.py --rscript /opt/R/4.4/bin/Rscript \
  -p proteins.fa --gene genes.txt -o out
```

Each run prints the resolved paths before doing any work:

```
External dependency check:
  Rscript ✔ (/opt/R/4.4/bin/Rscript)
  emapper ✔ (/opt/eggnog/emapper.py)
```

A missing or unusable required tool aborts the run with a non-zero status and a
message naming the flag to use. `Onestep` forwards its own overrides to the five
steps it runs, so a single set of flags is enough for a full analysis.

## 3. Temp and index hygiene

- `Target_finder` writes `<output>/temp` while scanning and deletes it at the end
  (previously the `rmtree` call was commented out, leaving CDS alignments behind).
  `--keep-temp` restores the old behaviour for debugging.
- `identification --clean` now deletes only `<prefix>/temp`. The genome index under
  `<output>/index` survives unless `--clean-index` is given; when the index came from
  `-d` it is never touched, and the skip is logged.

## 4. Restored files

- `src/miRNA_consistency.py` (544 lines) — three modes: rename against a reference
  basic-info (`-b`, no `-c`), unify non-conserved family names across groups
  (`-b -c`), and plain merge of several basic-info files (no `-b`).
- `scripts/dedup_score_pair.py` — keeps the highest-scoring row per unordered
  `(query, subject)` pair, as `docs/core-dataset.md` assumes.

## 5. Validation

- `./test.sh -t 12 -k` → `ALL TESTS PASSED` (identification, annotation, both
  target engines, determinism).
- Functional smoke suite, 33 checks: help text of the six analysis subcommands,
  `--step` parsing equivalence (`-h --step 1,2` byte-identical to `--step 1,2 -h`),
  both `miRNA_consistency` modes on a real basic-info file, `Stat`, `TFBS` on a BED,
  `Target_finder` with and without `--keep-temp` (outputs identical), the ssearch
  engine, early abort on a bad `--rscript` / `--emapper`, and `Differential_expression`
  on a negative-binomial 3v3 matrix (237 miRNAs, 33 up / 22 down, temp cleaned).
- Note: a first DE attempt used an 18-row synthetic matrix with near-constant
  dispersion; DESeq2 refuses to fit such data (`estimateDispersions` error). That is a
  property of the toy input, not of the code — the negative-binomial matrix above runs
  clean.

## 6. Upgrade notes

- No output format changed. `Target_finder` results are unchanged; only the temp
  directory behaviour differs.
- Scripts that relied on `identification --clean` deleting the index must now add
  `--clean-index`.

## 7. GitHub sync checklist

1. `git add -A && git commit` in the release tree; tag `v3.1.7c`.
2. Push the commit and the tag to `YangXZ-lab/mirdeep-p3`.
3. If a Docker image is rebuilt, tag it `3.1.7c-full` per `.github/workflows/docker-release.yml`.
4. `CHANGELOG.md` already carries the `[Unreleased] - 3.1.7c` section; move it under a
   dated heading when the tag is cut.
