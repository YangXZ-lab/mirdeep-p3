# miRDeep-P3 3.1.8b — release notes

_Released: 2026-11-18 · previous release: 3.1.8a · type: housekeeping (no algorithm or default changed)_

## 1. Summary

3.1.8b closes four of the improvement items in `docs/MIRDEEP-P3-ISSUES.md` (§3, O-xx):
the coordinate arithmetic is now one tested module, the combined `--step` help is generated
from the real parsers instead of being hand-maintained, the tree can be installed as a
normal console script, and the installed data index finally states which release it is.
identification/annotation results are byte-identical to 3.1.8a on the test dataset.

## 2. What changed

| Item | What | Where |
|---|---|---|
| O-05 | shared, unit-tested coordinate functions | `src/utils/coords.py` (new), `tests/test_coordinates.py` (new) |
| O-06 | combined `--step` help generated from parsers | `src/utils/help_display.py`, `mirdeep-p3`, `src/commands/analysis/__init__.py` |
| O-09 | standard console entry point (conservative) | `pyproject.toml`, `mirdeep_p3/__init__.py`, `mirdeep_p3/cli.py` (new) |
| O-10 | data-index release is reported, not enforced | `src/utils/data_index.py` (new), `tests/test_data_index.py` (new) |
| — | `--step` help/unit checks in the test suite | `test.sh` section 11, exit code 7 |
| — | `MIRDEEP_P3_HOME` honoured by `./mirdeep-p3` | `mirdeep-p3` |

## 3. Why the data index is reported but never enforced (O-10)

The software and the data index are released on two independent schedules, and the software
moves faster. That is exactly why a hard comparison against the software version would be
useless: `software 3.1.8b + index v3.1.4c` is the normal state, so the check would warn on
every run and be switched off within a week.

The real risk is different: the index silently ages while nobody notices, the ncRNA masking
starts to differ from the published results, and a run can no longer be reproduced. So 3.1.8b
compares **the installed index tag against `DATA_INDEX_RELEASE` in `data-index.env`** — both
sit on the data axis — and only reports:

- one line at start-up: `[index] data index: mirdeep-p3-v3.1.4c-full (matches data-index.env)`
- the full state is written to the batch manifest (`params.data_index`) and the installed tag
  to each per-sample manifest
- a changed tag invalidates `--resume`, because different masking means different reads
- a mismatch never aborts a run; only a **missing** index can do that

## 4. O-05 · coordinate arithmetic in one place

`src/utils/coords.py` now owns the two functions that had been copy-pasted inside
`src/primary-basic-info.py` and `src/mod_rm_redundant_meet_plant.py` (116 lines removed):

- `compute_star(pre_seq, mature_seq, pre_start, pre_end, strand, star_len=21)` returns a
  half-open `(start, end)` whose length always equals the returned sequence length — the exact
  invariant that 3.1.6c fixed for the reverse strand.
- `map_prediction(...)` turns relative hairpin coordinates into absolute genomic ones and
  raises `ValueError(negative or zero positions)`; the caller keeps its original warning text.

`tests/test_coordinates.py` (standard library only, runs under `python3 tests/test_coordinates.py`
or pytest) pins both strands × both mature placements against hand-computed numbers, the
star/mature non-overlap rule, and the swap when `beg > end`.

## 5. O-06 · combined help is now generated

Previously `--step 1,2 -h` used a 55-line hand-written help and every other combination fell
back to a placeholder that printed `(Full combined help will be implemented later.)`.

Now `mirdeep-p3` collects the real parsers of the selected steps and prints one table, grouping
options by which steps accept them (`Accepted by every selected step`, `Accepted by: 1, 2`),
marking options whose help text differs between steps, and listing the analysis subcommands.
`--step 1,2` additionally prints the combined run layout (output directories, auto-filled
`-i`, required `--prefix_miRNA` / `--species`, fallbacks from identification).

## 6. O-09 · console entry point (conservative)

`pip install -e .` now provides `mirdeep-p3` on `PATH` through `mirdeep_p3.cli:main`, which
finds the checkout and runs the existing entry script with `runpy`. Nothing was moved:
`src/`, `scripts/`, `data/` and the `./mirdeep-p3` script stay where they are, and only the
thin `mirdeep_p3` package is packaged. `MIRDEEP_P3_HOME` points the console script (and
`./mirdeep-p3`) at a tree that lives elsewhere.

## 7. O-10 · index status values

| State | Meaning | Effect |
|---|---|---|
| `missing` | `data/index` does not exist | identification/annotation cannot mask ncRNA |
| `unknown` | no `.index-version` and no `DATA_INDEX_RELEASE` | reported only |
| `unrecorded` | no `.index-version`, but `data-index.env` names a release | reported only |
| `match` | installed tag equals `DATA_INDEX_RELEASE` | nothing to do |
| `mismatch` | tags differ | reported only; the run is fine, but results may differ from a run made with the published index |

Write the tag with `echo mirdeep-p3-v3.1.4c-full > data/index/.index-version` (the conda and
Docker packages already do it; `test.sh` writes it after unpacking).

## 8. Validation

| Check | Result |
|---|---|
| `python3 tests/test_coordinates.py` | all checks passed |
| `python3 tests/test_data_index.py` | 5/5 passed |
| `--step 1,2` / `1,3` / `2,3` / `1,2,3` `-h` | all exit 0; `1,3` has no `--prefix_miRNA`, `2,3` does, `1,2,3` prints all three groups |
| `python3 -m mirdeep_p3.cli -h` | banner and help printed |
| identification smoke (`/nfs/users/zhaojw/Genburten/temp/idx-run`) | `[index] ... release unknown ...`, 12/12 stages, batch and sample manifests record `data_index` |
| `--resume` after the change | 12/12 stages skipped, new batch manifest with a second-resolution name |
| annotation smoke | `[index]` line printed, basic-info / cluster / count / exp all produced |
| `./test.sh -t 12 -k` | ALL TESTS PASSED (identification, annotation, Target_finder psrna + ssearch, Stat, DE, CLI/help/unit section) |
| identification outputs vs 3.1.8a | identical md5 for `mini_nr_predictions`, `.bed`, `mini_filter_P_prediction`, `mini.processed.fa`, `mini.total_reads` |

## 9. Upgrade notes

- No command-line option was removed or changed. `--step` combinations that used to print the
  placeholder now print real help.
- If you script against the generated help text, note that the layout changed.
- `mirdeep-p3` still runs without `pip install -e .`; the entry point is optional.
- An index that was already unpacked without `.index-version` will report `unrecorded` once
  and then stay quiet until you write the marker — nothing else changes.

## 10. GitHub sync checklist (when this release is pushed)

1. `git push origin main` from `/nfs/users/zhaojw/software/miRDeep-p3-3.1.7a` (no remote is
   configured yet; add `https://github.com/YangXZ-lab/mirdeep-p3.git` first).
2. Tag `v3.1.8b`; `data-index.env` needs no change (`mirdeep-p3-v3.1.4c-full` is still the
   published index).
3. Confirm `mirdeep-p3` is committed as mode `100755`.
4. Re-check that `docs/MIRDEEP-P3-INTERNALS.md` §3.2 / §4.4.1 / §11.1 / §12.3 / §12.4 / §13 and
   `docs/MIRDEEP-P3-ISSUES.md` §3 / §5 carry the 3.1.8b commit hash.
