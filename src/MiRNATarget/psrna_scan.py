#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
psrna_scan.py -- psRNATarget-compatible miRNA target scanner.

Pure-Python reimplementation of the psRNATarget (Dai et al., NAR 2018, V2) target
search, independent of ssearch36.  It is used by 'mirdeep-p3 analysis Target_finder
--engine psrna' and can also be run standalone.

Scoring rules (reverse-engineered and verified against archived psRNATarget output
on four plant datasets; see docs/PSRNA-ENGINE.md):

  * only positions 1..19 of the mature miRNA are scored (20 nt and longer tails
    do not contribute);
  * Watson-Crick pair 0.0; G:U wobble 0.5 (NOT up-weighted in the seed);
  * mismatch 1.5 inside the seed (positions 2..13), 1.0 outside;
  * single-base gap 2.0, plus 1.0 extra when the gap is on the miRNA side;
  * Expectation = sum of the above; sites with value > -e are discarded;
  * Inhibition = Translation when position 10 or 11 carries a real mismatch,
    otherwise Cleavage;
  * Multiplicity = number of non-overlapping sites of this miRNA on this target.

Output is the native psRNATarget 14-column table (tab separated, one header line).

Search strategy: all targets are concatenated into one string separated by 25 Ns,
encoded as a numpy uint8 array.  A gapless pass computes the penalty of every
position with a 19-mer sliding comparison, and a second pass anchors 8-9 nt exact
k-mers to recover single-base-indel candidates.  Each candidate is then scored
exactly (gapless, one deletion, one insertion) and the best non-overlapping sites
per (miRNA, target) pair are reported.
"""

import argparse
import os
import sys
from collections import Counter

try:
    import numpy as np
except ImportError:  # pragma: no cover
    sys.stderr.write("[psrna] numpy is required (conda: numpy>=1.26)\n")
    raise

ENC = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3, "N": 4}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "U": "A", "N": "N"}
SEP = "N" * 25
G = {}


def clean_name(name):
    """Strip trailing punctuation the way psRNATarget does (keeps duplicates)."""
    name = name.strip()
    while name and not name[-1].isalnum():
        name = name[:-1]
    return name


def read_fasta(path):
    """Return (names, sequences); sequences are upper-case DNA (U -> T)."""
    names, seqs, name, buf = [], [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    names.append(name)
                    seqs.append("".join(buf))
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.upper().replace("U", "T"))
    if name is not None:
        names.append(name)
        seqs.append("".join(buf))
    out_names, out_seqs = [], []
    for nm, sq in zip(names, seqs):
        nm = clean_name(nm)
        if not nm:
            continue
        out_names.append(nm)
        out_seqs.append(sq)
    return out_names, out_seqs


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGTUN", "TGCAAN"))


def is_gu(a, b):
    """G:U wobble between a miRNA base and the paired target-strand base."""
    return (a == "G" and b in "TU") or (a in "TU" and b == "G")


def score_pair(mi, tg, pos):
    """Penalty of one alignment column; pos is 1-based in miRNA orientation."""
    if pos > 19:
        return 0.0
    if mi == "-" or tg == "-":
        return 2.0
    if tg == COMP.get(mi):
        return 0.0
    if is_gu(mi, tg):
        return 0.5
    return 1.5 if 2 <= pos <= 13 else 1.0


def _align(mi, tgt):
    """Align miRNA mi against the target strand tgt (both 5'->3', '-' allowed)."""
    rev = tgt[::-1]
    cols = []
    i = 0
    pos = 0
    for ch in rev:
        mb = mi[i] if i < len(mi) else "-"
        i += 1
        if ch == "-":
            if mb != "-":
                pos += 1
            cols.append((mb, "-", pos))
        elif mb == "-":
            cols.append(("-", ch, pos))
        else:
            pos += 1
            cols.append((mb, ch, pos))
    while i < len(mi):
        mb = mi[i]
        i += 1
        if mb != "-":
            pos += 1
        cols.append((mb, "-", pos))
    return cols


def cost(cols):
    total = 0.0
    for mb, tb, pos in cols:
        value = score_pair(mb, tb, pos)
        if mb == "-" and tb != "-":
            value += 1.0
        total += value
    return total


def inhibition(cols):
    for mb, tb, pos in cols:
        if mb == "-" or tb == "-":
            continue
        if tb == COMP.get(mb) or is_gu(mb, tb):
            continue
        if pos in (10, 11):
            return "Translation"
    return "Cleavage"


def render(cols):
    """Return (miRNA fragment, alignment string, target fragment) for psRNATarget."""
    mf = "".join(c[0] for c in cols)
    tf = "".join(c[1] for c in cols)
    chars = []
    for mb, tb, pos in cols:
        if mb == "-" or tb == "-":
            chars.append(" ")
        elif tb == COMP.get(mb):
            chars.append(":")
        elif is_gu(mb, tb):
            chars.append(".")
        else:
            chars.append(" ")
    return (mf.replace("T", "U"),
            "".join(reversed(chars)),
            tf.replace("T", "U")[::-1])


def best_site(mi, target, offset, cutoff):
    """Score gapless / one-gap alignments at one position; return best or None."""
    length = len(mi)
    best = None
    trials = []
    if offset + length <= len(target):
        trials.append(("plain", target[offset:offset + length], length, None))
    if G.get("allow_gap", False) and length >= 2 and offset + length - 1 <= len(target):
        site = target[offset:offset + length - 1]
        for k in range(length):
            trials.append(("del", site[:k] + "-" + site[k:], length - 1, k))
    if G.get("allow_gap", False) and offset + length + 1 <= len(target):
        site = target[offset:offset + length + 1]
        for k in range(length + 1):
            trials.append(("ins", site, length + 1, k))
    for kind, site, span, k in trials:
        if kind == "ins":
            cols = _align(mi[:k] + "-" + mi[k:], site)
        else:
            cols = _align(mi, site)
        penalty = cost(cols)
        if penalty <= cutoff + 1e-9:
            if kind == "plain":
                return (penalty, cols, span)
            if best is None or penalty < best[0] - 1e-9:
                best = (penalty, cols, span)
    return best


def scan_one(idx):
    """Scan every target for one miRNA; returns a list of site tuples."""
    code_arr = G["S"]
    big = G["big"]
    starts = G["starts"]
    tseqs = G["tseqs"]
    tnames = G["tnames"]
    mname, mseq = G["mirnas"][idx]
    cutoff = G["cutoff"]
    length = len(mseq)
    query = revcomp(mseq)
    klen = min(19, length)
    j0 = length - klen
    core = query[j0:]
    n = code_arr.size - klen + 1
    penalty = np.zeros(n, dtype=np.uint8)
    buf = np.empty(n, dtype=np.uint8)
    for k in range(klen):
        base = core[k]
        pos = klen - k
        weight = 3 if 2 <= pos <= 13 else 2
        lut = np.full(5, weight, dtype=np.uint8)
        lut[ENC[base]] = 0
        if base == "A":
            lut[ENC["G"]] = 1
        elif base == "C":
            lut[ENC["T"]] = 1
        np.take(lut, code_arr[k:k + n], out=buf, mode="clip")
        np.add(penalty, buf, out=penalty)
    cand = [int(i) - j0 for i in np.flatnonzero(penalty <= int(round(cutoff * 2)))]
    if klen >= 19:
        anchors = [(s, 9) for s in (0, 10)] + [(s, 8) for s in (4, 8, 11)]
        for s0, alen in anchors:
            seed = core[s0:s0 + alen]
            hit = big.find(seed)
            while hit != -1:
                base = hit - (j0 + s0)
                cand.append(base)
                cand.append(base - 1)
                cand.append(base + 1)
                hit = big.find(seed, hit + 1)
    rows = []
    seen = set()
    for p in cand:
        if p in seen or p < 0:
            continue
        seen.add(p)
        si = int(np.searchsorted(starts, p, side="right")) - 1
        if si < 0:
            continue
        offset = p - int(starts[si])
        target = tseqs[si]
        if offset < 0 or offset >= len(target):
            continue
        res = best_site(mseq, target, offset, cutoff)
        if res is None:
            continue
        penalty_value, cols, span = res
        mf, aln, tf = render(cols)
        rows.append((idx, si, offset, penalty_value, mname, tnames[si],
                     offset + 1, offset + span, mf, aln, tf, inhibition(cols)))
    rows.sort(key=lambda r: (r[1], r[3], r[2]))
    kept = []
    used = {}
    for row in rows:
        lo, hi = row[6], row[7]
        overlap = False
        for (a, b) in used.get(row[1], ()):
            if lo <= b and a <= hi:
                overlap = True
                break
        if overlap:
            continue
        used.setdefault(row[1], []).append((lo, hi))
        kept.append(row)
    kept.sort(key=lambda r: r[2])
    return kept


def add_arguments(parser):
    parser.add_argument("-m", "--mirna", required=True,
                        help="Mature miRNA FASTA file")
    parser.add_argument("-t", "--targets", required=True,
                        help="Target (CDS / transcript) FASTA file")
    parser.add_argument("-o", "--output", required=True,
                        help="Output table (native psRNATarget 14-column TSV)")
    parser.add_argument("-e", "--expect", type=float, default=2.5,
                        help="Expectation cutoff (default: 2.5)")
    parser.add_argument("-p", "--threads", type=int, default=1,
                        help="Number of worker processes (default: 1)")
    parser.add_argument("--allow-gap", action=argparse.BooleanOptionalAction, default=True,
                        help="Report single-base-indel sites too (default: on; "
                             "--no-allow-gap keeps only indel-free sites)")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="psRNATarget-compatible miRNA target scanner (no ssearch36 needed).")
    add_arguments(parser)
    args = parser.parse_args(argv)

    tnames, tseqs = read_fasta(args.targets)
    mnames, mseqs = read_fasta(args.mirna)
    if not tseqs:
        sys.exit("[psrna] no target sequence found in " + args.targets)
    if not mseqs:
        sys.exit("[psrna] no miRNA sequence found in " + args.mirna)
    sys.stderr.write("[psrna] %d targets / %d bases, %d miRNAs, expect <= %.1f\n"
                     % (len(tseqs), sum(len(s) for s in tseqs), len(mnames), args.expect))

    G["big"] = SEP.join(tseqs)
    starts, offset = [], 0
    for seq in tseqs:
        starts.append(offset)
        offset += len(seq) + len(SEP)
    G["starts"] = np.array(starts, dtype=np.int64)
    raw = np.frombuffer(G["big"].encode("ascii"), dtype=np.uint8)
    code_table = np.full(256, 4, dtype=np.uint8)
    for base, value in ENC.items():
        code_table[ord(base)] = value
    G["S"] = code_table[raw]
    G["tseqs"] = tseqs
    G["tnames"] = tnames
    G["mirnas"] = list(zip(mnames, mseqs))
    G["cutoff"] = args.expect
    G["allow_gap"] = bool(args.allow_gap)

    rows = []
    if args.threads > 1:
        import multiprocessing as mp
        workers = min(args.threads, len(mnames))
        with mp.Pool(workers) as pool:
            for i, res in enumerate(pool.imap_unordered(scan_one, range(len(mnames))), 1):
                rows.extend(res)
                if i % 20 == 0:
                    sys.stderr.write("[psrna] %d/%d miRNAs, %d sites\n"
                                     % (i, len(mnames), len(rows)))
    else:
        for i in range(len(mnames)):
            rows.extend(scan_one(i))
            if (i + 1) % 20 == 0:
                sys.stderr.write("[psrna] %d/%d miRNAs, %d sites\n"
                                 % (i + 1, len(mnames), len(rows)))

    rows.sort(key=lambda r: (r[4], r[5], r[6]))
    mult = Counter((r[0], r[1]) for r in rows)
    out_dir = os.path.dirname(os.path.abspath(args.output))
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    with open(args.output, "w") as out:
        out.write("miRNA_Acc.\tTarget_Acc.\tExpectation\tUPE$\tmiRNA_start\tmiRNA_end\t"
                  "Target_start\tTarget_end\tmiRNA_aligned_fragment\talignment\t"
                  "Target_aligned_fragment\tInhibition\tTarget_Desc.\tMultiplicity\n")
        for row in rows:
            idx, si, _start, penalty_value, mname, tname, ts, te, mf, aln, tf, inhib = row
            out.write("%s\t%s\t%.1f\t-1.0\t1\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t\t%d\n"
                      % (mname, tname, penalty_value, len(mf.replace("-", "")),
                         ts, te, mf, aln, tf, inhib, mult[(idx, si)]))
    sys.stderr.write("[psrna] total sites: %d -> %s\n" % (len(rows), args.output))


if __name__ == "__main__":
    main()
