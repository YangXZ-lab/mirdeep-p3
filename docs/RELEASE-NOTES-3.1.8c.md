# miRDeep-P3 3.1.8c — release notes

_Released: 2026-11-18 · previous release: 3.1.8b · type: correctness (opt-in correction)_

## 1. Summary

3.1.8c fixes a silent 1-nt error in the mature 5' boundary. The mature interval
reaches the pipeline from the prediction and is then adjusted per arm by
`utils/coords.py`; that adjustment can leave it one base longer or shorter at the
5' end than the reads that support the locus. Nothing fails when it happens, and
the expression matrix absorbs it as a 1-nt offset instead of reporting it.

The release adds an audit that measures the offset against the reads on every run,
and an opt-in correction (`--fix-mature-boundary`) that rewrites the table only
when the reads disagree with it. Defaults leave every product byte-identical to
3.1.8b.

## 2. What changed

| Item | What | Where |
|---|---|---|
| — | coordinate arithmetic for the boundary problem, stdlib only | `src/utils/mature_boundary.py` (new) |
| — | CLI: reads `*.mature.aln`, writes `<prefix>-boundary-audit.tsv` | `src/audit_mature_boundary.py` (new) |
| — | unit tests, incl. a real regression case | `tests/test_mature_boundary.py` (new) |
| — | opt-in correction | `annotation --fix-mature-boundary` |
| — | evidence thresholds (0.60 / 10.0) | `annotation --boundary-min-fraction`, `--boundary-min-weight` |
| — | part 2 order: align → audit (+correct) → cluster → expression | `src/commands/annotation.py` |
| — | version bump and release notes | `src/version.py`, `CHANGELOG.md`, `README.md`, `docs/DATA-INDEX.md`, `test.sh` |

## 3. Why the audit reports instead of correcting

Measured on four annotated datasets — Arabidopsis 135, Gnetum 513, rice 446,
Scaevola 344 records; 1420 loci with usable coverage — 47 loci (3.3%) have a
dominant read offset that is not canonical. Three are confirmed annotation errors
against miRBase: `Nip-MIR166a`, `Nip-MIR166b` (the annotation keeps one extra base
at the 5' end; the reads are the miRBase mature) and `Gmo-MIR159a`. The first two
are the most abundant members of that family at 1.6M counts each, so the error is
not numerically small.

What makes the diagnosis awkward: 24 of the 47 are real 5'-truncated molecules,
not annotation errors. Their annotated mature already sits on the arm boundary, so
there is no upstream room for the annotation to have drifted. The audit therefore
reports rather than assumes, and the correction is opt-in.

## 4. trim or shift — the distinction that matters

An annotated mature that disagrees with its dominant read can be wrong in two
different ways, and they need different edits:

| Locus | Annotated mature | Dominant read | Correct action |
|---|---|---|---|
| `Nip-MIR166a` | `TTCGGACCAGGCTTCATTCCC` (21) | `TCGGACCAGGCTTCATTCCCC` (21) | **shift** — drop one base at 5', add one at 3' |
| `Gmo-MIR159a` | `GTTGGATTGAAGGGAGCTCCA` (21) | `TTGGATTGAAGGGAGCTCCA` (20) | **trim** — drop one base at 5' |

Trimming `MIR166a` would produce a 20-nt mature, further from miRBase than the
original. The implementation therefore builds both candidates and accepts only the
one that reproduces the dominant read sequence exactly; if neither does, the row
is left alone.

## 5. Two windows, and they are not interchangeable

A 22-column basic-info carries two different sequence windows and mixing them
silently disables the containment check:

- columns 7-9: the precursor ±20 nt window
- columns 11-13: the real hairpin (mature + loop + star)

Candidates are built in the wide window — it is wide enough that a mature shifted
by one base still lands inside it, which the real hairpin is not. The
arm-placement label and the containment check use the real hairpin. Reading only
the wide window made the arm label a 20-nt-too-large yardstick; that bug was found
and fixed during this release, and both column layouts (22 and 24) now agree on
every locus.

## 6. Guards

- **Only ±1 nt corrections are applied.** A locus off by two or more is reported
  as `shift_beyond_limit` and left alone — at that distance the dominant offset is
  more likely 5'-degradation than a boundary slip.
- **No automatic upstream extension.** Reads starting before the annotated mature
  are reported as `shift_upstream`; rewriting the boundary upstream on that
  evidence is a much stronger claim than removing a base the reads do not support.
  `--allow-extend` exists for the case where it is wanted.
- **A correction touches the mature triplet and nothing else.**

## 7. Known limitation

A correction can push the mature one base past the annotated precursor end, which
is what happens when the hairpin was excised flush with the old (wrong) 3' end:
`Nip-MIR166a` and `Nip-MIR166b` on rice are exactly that. They stay corrected,
because the reads are the evidence, and the run log names them:

```
[boundary]   corrected mature overhangs the annotated precursor: 2 (Nip-MIR166a, Nip-MIR166b)
```

Extending the precursor would also mean re-folding the hairpin and rewriting the
stored structure, so this is reported rather than edited silently. 3.1.8d closes
this gap.

## 8. Validation

| Check | Result |
|---|---|
| `python3 tests/test_mature_boundary.py` | all checks passed (local and on the server) |
| `./test.sh -k -t 8` | ALL TESTS PASSED |
| default products vs 3.1.8b | basic-info `f9223f3c...`, cluster `db799e4e...` — identical |
| rice dataset, audit | 446 loci: 405 canonical / 18 low coverage / 23 flagged (2 shift, 1 upstream, 12 beyond limit, 7 ambiguous, 1 out of range) |
| rice dataset, correction | 2 rows rewritten (`Nip-MIR166a`, `Nip-MIR166b`); corrected sequence equals miRBase `osa-miR166a-3p` |
| corrected row shape | 446 rows × 22 columns = 9812 cells; 6 cells change, in 2 rows, columns 16/17/18 only |
| clustering | 0-line difference between original and corrected table |
| expression reference | 446 entries, same order, 2 differ, length unchanged |
| expression matrices | 2 of 446 rows change value; the other 444 are identical value-for-value |
| fixture end-to-end (reads shifted by 1 nt) | 19/19 flagged → 19 corrected → post-fix audit 19/19 canonical |

## 9. Upgrade notes

- No command-line option was removed or changed. The audit runs every time and
  costs seconds; it prints one summary line per run.
- `--fix-mature-boundary` is off by default. Turning it on changes at most a
  handful of rows per dataset, and only the mature triplet.
- The annotation part 2 order changed (clustering now runs after the audit). Any
  script that called the clustering step expecting the pre-audit table now sees
  the audited/corrected one — which is the intent.

## 10. Verification record

The full measurement log, including the four-species scope check and the
regression fixture, is in `docs/verification/TEST-3.1.8c.md`.
