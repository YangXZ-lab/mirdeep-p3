# miRDeep-P3 3.1.8d — release notes

_Released: 2026-11-18 · previous release: 3.1.8c · type: correctness (completes the 3.1.8c correction)_

## 1. Summary

3.1.8d closes the loop that 3.1.8c left open. A corrected mature can end up one
base outside its own precursor, and until now that was reported but not fixed.
The new step widens the precursor on the short side, reads the missing bases from
the genome, re-folds the widened precursor and its window with RNAfold, and
recomputes the star.

It runs only when `--fix-mature-boundary` actually rewrote the table, and it never
touches the mature sequence — so the bowtie reference and the alignments stay
valid and no third round of alignment is needed.

## 2. What changed

| Item | What | Where |
|---|---|---|
| — | random access to a genome FASTA without writing `.fai` | `src/utils/genome_fasta.py` (new) |
| — | widening decision, transcript-offset check, guard rails, star recomputation | `src/utils/precursor_repair.py` (new) |
| — | CLI, also usable standalone | `src/repair_precursor.py` (new) |
| — | unit tests | `tests/test_precursor_repair.py` (new) |
| — | repair step in part 2, after correction and before clustering | `src/commands/annotation.py` |
| — | basic-info interval documented as 0-based (3.1.8c said 1-based) | `src/utils/mature_boundary.py`, `src/audit_mature_boundary.py` |
| — | version bump, release notes, test suite hook | `src/version.py`, `CHANGELOG.md`, `README.md`, `docs/DATA-INDEX.md`, `test.sh` |

## 3. The problem in one picture

The precursor is excised flush with the mature it was built around. Move the
mature one base downstream and its 3' end lands one base outside the hairpin:

```
Nip-MIR166a   mature   3487877..3487898
              precursor 3487799..3487897     short by 1 nt at the high end
```

Nothing crashes, but the mature stops being a substring of the precursor, so the
stored dot-bracket describes a window that no longer contains the mature it
annotates, and the star is placed relative to a frame that is now off by one.

Measured on the rice dataset: 0 of 446 loci overhang before the correction, 2
after (`Nip-MIR166a`, `Nip-MIR166b`), both by exactly 1 nt. Of 1420 loci across
four species, only those two are affected.

## 4. The genome check settles the coordinate convention

The repair is the first code in the tree that slices a genome with basic-info
coordinates, so the convention had to be established rather than assumed. 3.1.8c
called the interval "1-based start, exclusive end". Testing all three readings
against IRGSP-1.0 on the 446 stored records:

| Convention | Precursor reproduces | Mature reproduces |
|---|---|---|
| 1-based start + exclusive end | 0 / 446 | 0 / 446 |
| **0-based start + exclusive end** | **446 / 446** | **446 / 446** |
| 1-based, both ends inclusive | 0 / 446 | 0 / 446 |

Only the 0-based reading reproduces the stored sequences, and the Arabidopsis test
dataset agrees (19/19). Nothing computed wrong before — only differences of two
coordinates are used — but the documentation was wrong and is corrected here.

## 5. Widening is verified to have only added

The repair does not trust the coordinates blindly. Before writing, it checks that
the old precursor is still present in the widened sequence at the offset the
coordinates imply. If it is not, the genome is a different assembly and the row is
left alone. The guard rails:

| Guard | Behaviour |
|---|---|
| `--max-extension` (default 5 nt per side) | refuse and report if more is needed |
| `sequence_mismatch` | old precursor not found at the implied offset → leave the row alone |
| `fold_failed` | RNAfold returns nothing → leave the row alone |
| `out_of_genome` | widened interval leaves a chromosome → leave the row alone |
| `star_unresolved` | widened but `compute_star()` returns nothing → keep the original star and report |

## 6. The star does not move, and that is not a coincidence

Both rice loci overhang at the **transcript 3' end**, while the star sits on the
opposite arm. Widening the 3' end therefore cannot move it, and the recomputed
star matches the original for every widened locus — including the 14 loci of the
end-to-end fixture. A path that would move it exists (`--allow-extend` on the 5'
side) and is implemented, but no measured locus triggers it.

## 7. Validation

| Check | Result |
|---|---|
| `python3 tests/test_precursor_repair.py` | all checks passed (local and on the server) |
| `./test.sh -k -t 8` | ALL TESTS PASSED |
| default products vs 3.1.8b/3.1.8c | basic-info `f9223f3c...`, cluster `db799e4e...` — identical; no repair report produced |
| rice dataset, corrected → repaired | 2 rows change, columns 6-13 only; the mature triplet changes in 0 of 446 rows |
| rice dataset, final table | mature overhanging the precursor: 0; mature not a substring: 0; structure length mismatch: 0 |
| fixture end-to-end (reads shifted by 1 nt, 19 loci) | 14/19 overhang → all widened and re-folded → post-repair audit 18 canonical + 1 ambiguous |
| star movement | 0 of 14 widened loci move their star |

## 8. Notes on `src/utils/genome_fasta.py`

`samtools faidx` and `bedtools getfasta` both want to write an index next to the
FASTA, which is not always allowed; pyfaidx would be a new dependency. This module
builds its own byte-offset index, writes it beside the file when that is writable
and falls back to a temp directory when it is not. A 360 MB genome costs one scan
(about 2 s) and a few hundred bytes per query.

## 9. Upgrade notes

- No command-line option was removed or changed. `--fix-mature-boundary` now also
  repairs precursors, and writes a second report, `<prefix>-precursor-repair.tsv`.
- The pre-repair table is kept as
  `<prefix>-basic-info.before-precursor-repair`.
- Anyone slicing a genome with basic-info coordinates must use the 0-based
  reading. The documentation in `mature_boundary.py` and
  `audit_mature_boundary.py` now says so.

## 10. Verification record

The measurement log, including the four-species scope check and the shifted-reads
fixture, is in `docs/verification/TEST-3.1.8d.md`.
