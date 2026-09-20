# The psRNATarget-compatible engine (--engine psrna)

> Applies to: miRDeep-P3 3.1.7b and later ｜ Code: src/MiRNATarget/psrna_scan.py
> Entry point: mirdeep-p3 analysis Target_finder (psrna is the default engine since 3.1.7b)
> Related: docs/RELEASE-NOTES-3.1.7b.md (what changed and why)

## 1. Why a second engine

Target_finder used to have exactly one engine: ssearch36 plus the jtremblay/mirnatarget
parser (scripts/parse_ssearch.py + src/MiRNATarget/parse_mirna_targets.py). That pipeline
scores with psRNATarget rules, but it only ever scores the HSPs ssearch36 reports, and
ssearch36 -b 200 caps how many HSPs a single miRNA can return. The result drifts away
from what psRNATarget actually produced.

Measured on Arctium_lappa (562 archived psRNATarget miRNA-target pairs):

| engine | pairs returned | shared with psRNATarget | missing | extra |
|---|---|---|---|---|
| ssearch (historical) | 605 | 352 (62.6%) | 210 | 253 |
| psrna (new) | 561 | 559 (99.5%) | 3 | 2 |

The new engine implements the same scoring rules but searches every position of every
target, so it neither truncates nor depends on an external aligner. Since 3.1.7b it is
the default engine, with `--allow-gap` on: with those defaults Arctium_lappa yields 579
rows (missing 1, extra 17, 99.82% identical with the alignment column excluded).

## 2. Scoring rules implemented

These were reverse-engineered from archived psRNATarget output and then verified: feeding
the archived alignment strings back through the scorer reproduces the archived
Expectation for 563/563 and 1723/1724 rows on two datasets, and a grid search over 1920
parameter combinations picked exactly this set as the unique optimum.

1. Only positions 1..19 of the mature miRNA are scored; bases 20 and beyond add nothing.
   A 20, 21 or 22 nt mature sequence therefore scores identically.
2. Per-position penalty:
   - Watson-Crick pair: 0.0
   - G:U wobble: 0.5 (NOT multiplied by the seed factor)
   - mismatch: 1.5 at positions 2..13 (seed), 1.0 elsewhere
   - single-base gap: 2.0 on either side, plus 1.0 when the gap sits on the miRNA side
3. Expectation = sum of the penalties. Sites above -e/--evalue (default 2.5) are dropped.
4. Inhibition: Translation when position 10 or 11 carries a real mismatch, else Cleavage.
5. Multiplicity: number of non-overlapping sites of that miRNA on that target. Sites of
   the same miRNA are forced to be non-overlapping; different miRNAs may share a site.
6. UPE is reported as -1.0 (accessibility is off, matching the archived output).
7. The alignment string is printed 3'->5' with respect to the mature miRNA, i.e. in the
   target's 5'->3' orientation: ':' pair, '.' G:U, ' ' mismatch or gap.

## 3. Search strategy

* All targets are concatenated into one string separated by 25 Ns and encoded as a numpy
  uint8 array. Every position of the concatenated string is then addressed at once.
* Gapless pass: the 19-mer of the miRNA (reverse complemented) is slid across the encoded
  array and the penalty of each position is accumulated column by column. Internally the
  weights are doubled (3 / 2 / 1) so the arithmetic stays integral; candidates are the
  positions with penalty <= 2 x cutoff.
* Indel pass (always on for anchoring, reporting gated by --allow-gap, on by default):
  8-9 nt exact
  k-mer anchors at several offsets locate positions where the target carries one extra or
  one missing base.
* Each candidate is scored exactly with three alignments (gapless, one deletion, one
  insertion) and the best non-overlapping sites per (miRNA, target) are kept.
* Work is distributed over processes with multiprocessing.Pool(-p/--threads); each worker
  gets one miRNA.

## 4. Validation against archived psRNATarget output

Four independent plant datasets; "identical, ignoring alignment" compares every field
except the alignment string, which psRNATarget renders only over the HSP extent.

| dataset | miRNAs | targets | psRNATarget rows | psrna rows | missing | extra | identical (ignoring aln) |
|---|---|---|---|---|---|---|---|
| Arctium_lappa | 60 | 44,686 | 563 | 562 | 3 | 2 | 560 / 563 = 99.47% |
| Helianthus_annuus | 143 | 47,137 | 1,724 | 1,729 | 3 | 8 | 1,716 / 1,724 = 99.54% |
| Lactuca_sativa | 211 | 47,803 | 2,531 | 2,588 | 7 | 64 | 2,522 / 2,531 = 99.64% |
| Bidens_alba | 511 | 165,431 | 12,296 | 12,776 | 151 | 631 | 12,105 / 12,296 = 98.45% |

Runtime: Arctium seconds, Lactuca about 15 min at -t 16, Bidens about 40 min at -t 32
(511 miRNAs x 165,431 CDS = 181 Mb).

### Where the remaining differences come from

* extra rows: psRNATarget only scores the HSPs ssearch reported, so a site that is optimal
  under the psRNATarget penalty but not an ssearch HSP never appears in its output. Proof:
  on Lactuca, Lsa-MIRN10c x XP_023768828.1 at 68..89 is an extra here, and ssearch36 on
  that single pair only reports HSPs at 242-263 / 1956-1975 / 125-144 / 2451-2471.
  psRNATarget also seems to cap rows per miRNA: on Bidens the four identical MIRN361
  members each stop at exactly 297 rows, while the exhaustive search finds 375.
* missing rows: mostly single-base indels (75 of 151 on Bidens), which are excluded unless
  --allow-gap is given; the rest are the same site reported one base off.
* alignment column: psRNATarget renders only the HSP extent, psrna renders all 19 nt. The
  score is unaffected; the tails differ.

## 5. CLI

Standalone:

    src/MiRNATarget/psrna_scan.py -m mature.fa -t cds.fa -o out.tsv -e 2.5 -p 32 [--no-allow-gap]

Through miRDeep-P3:

    mirdeep-p3 analysis Target_finder -i mature.fa -c cds.fa -o out -t 32   # psrna is the default

Options that only apply to this engine: --allow-gap / --no-allow-gap, --psrna-detail.
Options that are ignored by this engine: --GUs (psRNATarget V2 has no G:U cap),
--max-hsp (ssearch36 specific). Both print a note or are rejected with an error.

## 6. Known limits

* Indel-containing sites are reported by default (`--allow-gap`, since 3.1.7b);
  `--no-allow-gap` excludes them. Keeping them costs 0.0-4.3% extra rows across the four
  validation datasets and recovers pairs that are otherwise missing (see section 1).
* Duplicate names in the mature FASTA are kept as separate scoring jobs, but the
  mature-length lookup used for the six-column table keeps the last occurrence.
* The engine needs numpy (already in mirdp3_environment.yml, numpy>=1.26). Nothing else.
