# miRDeep-P3 3.1.8e — release notes

_Released: 2026-11-18 · previous published release: 3.1.6c · type: correctness (drops loci) + accumulated features_

## 1. Summary

3.1.8e adds the check that was missing between "the mature is placed on the
precursor" and "the precursor is a hairpin": the two arms have to leave room for a
loop. When the precursor is a window cut tight around the mature, they do not, and
a locus can end up with a star that overlaps its own mature.

The filter runs inside annotation part 1, right after basic-info is built and
before any alignment, so no reads are aligned, no hairpin is folded and no matrix
is counted for a locus that will be discarded.

**This is the first release in the 3.1.8 line that changes the row count of
basic-info by default.** On the rice dataset that is 446 → 432.

## 2. Release history — why this release is 3.1.8e, not 3.1.6d

The repository's previous published release was **3.1.6c** (2026-09-18). Six
releases have been cut since then and are all present in this tree, in order:

| Version | Date | What it changed | Default output affected? |
|---|---|---|---|
| 3.1.6d | — | `-fa` zero-filling, `--max-5p-trunc` / `--max-3p-trunc`, `Onestep --min-expr`. Never published on its own; folded into 3.1.8b. | no |
| 3.1.7b | 2026-09-20 | psRNATarget-compatible engine becomes the default for `Target_finder`; test fixtures leave version control | target prediction |
| 3.1.7c | 2026-09-26 | unified external tool resolution, temp cleanup, restored consistency script | no |
| 3.1.8a | 2026-10-29 | identification becomes 12 declarative stages; `--resume`; run manifests | no |
| 3.1.8b | 2026-10-29 | coordinate helpers in one tested module, generated `--step` help, console entry point, data-index reporting | no |
| 3.1.8c | 2026-11-18 | audit and opt-in correction of the mature 5' boundary | only with `--fix-mature-boundary` |
| 3.1.8d | 2026-11-18 | repair precursors that the corrected mature overhangs | only with `--fix-mature-boundary` |
| **3.1.8e** | **2026-11-18** | **arm-gap filter** | **yes — row count** |

Each step has its own release note in `docs/` and its own commit, so the
intermediate state can be inspected without guessing: `git log` walks
3.1.6c → 3.1.7b → 3.1.7c → 3.1.8a → 3.1.8b → 3.1.8c → 3.1.8d → 3.1.8e and each
commit's tree is the tree that was tested for that version.

Only 3.1.8e is published as an artefact (conda package, Docker image, git tag).
The intermediate versions were development checkpoints and were never shipped.

## 3. What changed

| Item | What | Where |
|---|---|---|
| — | the predicate: net arm gap, structure presence, coordinate sanity | `src/utils/structure_filter.py` (new) |
| — | CLI, usable standalone against any 22-column basic-info | `src/filter_structure.py` (new) |
| — | unit tests | `tests/test_structure_filter.py` (new) |
| — | filter runs at the end of part 1 | `src/commands/annotation.py` |
| — | threshold, default 10 nt | `annotation --min-arm-gap` |
| — | report, one row per input locus | `<prefix>-structure-filter.tsv` |
| — | version bump, release notes, test suite hook | `src/version.py`, `CHANGELOG.md`, `README.md`, `docs/DATA-INDEX.md`, `test.sh` |

## 4. The problem

`compute_star()` takes `star_len` bases from the opposite end of the precursor and
only verifies that the mature sits at *one* end. It never checks that the two arms
clear each other. When the precursor is a window cut tight around the mature, they
do not:

```
Nip-MIR12457a   precursor 3742249..3742289   40 nt
                mature    3742249..3742269   20 nt
                star      3742268..3742289   21 nt     20 + 21 = 41 > 40
```

The star starts inside the mature, and the stored dot-bracket describes a hairpin
that neither arm matches. Nothing downstream fails.

Measured on the rice dataset: 1 locus of 446 overlaps. Relaxing to "a loop of at
least 10 nt" — a real plant hairpin loop is 15-100+ nt — catches 14 (3.1%), and all
14 share the same signature: a precursor of 40-50 nt that is exactly
`mature + star + loop` (e.g. `MIR12389a`: 50 = 21 + 21 + 8), far below the p5 of
the dataset (64 nt).

## 5. These are not repairable by widening the precursor

3.1.8d added the ability to widen a precursor, so widening was the obvious
response. It does not work here. Every one of the 14 loci has its mature in the
middle of the ±20 nt window rather than on an arm, and `compute_star()` requires
the mature to sit at one end. Widening to the whole window makes it return `None`:
there is no hairpin with this mature as an arm. Dropping the locus is the only
correct action.

## 6. What the predicate actually tests

It tests whether the precursor has room for two arms and a loop. It does **not**
test whether the star is expressed — stars are usually inferred rather than
observed. `--min-arm-gap 0` reduces it to "do not overlap".

| Status | Condition |
|---|---|
| `arm_overlap` | gap < 0 |
| `short_arm_gap` | 0 ≤ gap < `--min-arm-gap` |
| `missing_precursor_structure` / `missing_stem_structure` | structure absent, or its length disagrees with the sequence |
| `missing_star` | star sequence is empty |
| `bad_coordinates` | non-integer, empty or inverted interval |

The filter only drops rows. It never edits one, so kept rows stay byte-identical.

## 7. Validation

| Check | Result |
|---|---|
| `python3 tests/test_structure_filter.py` | all checks passed (local and on the server) |
| `./test.sh -k -t 8` | ALL TESTS PASSED |
| bundled mini dataset | 19/19 pass, 0 dropped; basic-info `f9223f3c...` unchanged, so nothing is over-filtered |
| rice dataset | 446 → **432** (14 dropped: 1 overlap + 13 short gap); the 432 kept rows are byte-identical to the input; the smallest surviving gap is 13 nt |
| Scaevola dataset, full run (15 libraries, 1.3 GB genome) | 344 → **343** (1 dropped, `Sta-MIR3711a`: precursor 48 = 21 + 21 + 6) |
| execution order | log line 162 filter → line 169 the bowtie reference is built → line 245 boundary audit: dropped loci are never aligned |
| per-locus quality, 343 loci × 5 checks | coordinate length vs sequence length, structure length, genome reproduction of mature/star/precursor, containment, gap ≥ 10 — all pass |
| gap distribution (343 loci) | min 16 / median 82 / max 235 |
| downstream | `Stat --rnaplot` produced 343 structure images; DE (leaf vs root) retained 152, 36 up / 36 down |
| row alignment | basic-info 343 = `mature.count` 343 = `mature.exp` 343, same names in the same order |

## 8. The star coordinate bug that was still hiding

The Scaevola dataset had a previous annotation produced by an older 3.1.6c-line
build. Matching the two runs by genomic locus (not by name — the MIRN families are
renumbered between runs):

| | old result | 3.1.8e |
|---|---:|---:|
| loci | 344 | 343 |
| loci in common | — | 343 |
| `star_end − star_start` equals the star sequence length | **0 / 306** | **306 / 306** |
| star sequence / precursor sequence / precursor structure | — | **343/343 identical** |
| loci whose star overlaps the mature | **170** | **0** |

The old file stored a 20-nt interval for a 21-nt star, which is why 170 loci looked
like they had overlapping arms. Apart from that, the two runs agree byte for byte —
the coordinate fix changed the coordinates it was meant to change and nothing else.

## 9. Compatibility

**basic-info can now have fewer rows than `-all-mod_fp_prediction` entries.** Any
downstream script that assumed the two counts match must read
`<prefix>-structure-filter.tsv`, which lists every input locus with its status and
the reason it was dropped.

The pre-filter table is kept as
`temp/<prefix>-basic-info.before-structure-filter`, so nothing is lost.

## 10. Upgrade notes

- `--min-arm-gap` defaults to 10 nt. Use `0` to keep every non-overlapping locus.
- The filter is not optional and cannot be switched off entirely; it is the
  structural validity check the pipeline was missing. Its strictness is the
  parameter.
- Every other default is unchanged from 3.1.8d: `--fix-mature-boundary` is still
  off, and with no failing locus the output equals 3.1.8d.

## 11. Verification records

- `docs/verification/TEST-3.1.8e.md` — unit, end-to-end and rice-dataset results
- `docs/verification/TEST-FULL-3.1.8e.md` — the Scaevola full run: 15 libraries,
  1.3 GB genome, annotation + Stat + differential expression
- `docs/verification/TEST-3.1.8c.md`, `TEST-3.1.8d.md` — the two preceding steps

## 12. How this release is published

Three artefacts, all from the same tree:

| Channel | Identifier | Note |
|---|---|---|
| git tag | `v3.1.8e` | source only; does not trigger the image build |
| git tag + Release | `v3.1.8e-full` | triggers `docker-release.yml`; carries the full source archive |
| conda | `jaguares/mirdeep-p3=3.1.8e` | built from `conda.recipe/`, with `data/index` bundled |

The separate `v3.1.8e` / `v3.1.8e-full` pair exists because
`docker-release.yml` filters on `v*-full`, so a plain version tag cannot start a
40-minute image build by accident.

`data-index.env` still names `mirdeep-p3-v3.1.4c-full`. The data index is released
on its own schedule and this software release does not change it.
