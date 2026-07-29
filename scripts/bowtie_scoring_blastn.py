#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
Accurate scoring of miRNA alignments from BLASTn‑short output.

Rules (reference‑centric):
  - Left overhang region: difference where one side is '-' => 5' overhang (1.0)
                          both sides have bases (mismatch) => internal (1.0)
  - Internal aligned region: any difference => 1.0
  - Right overhang region: difference where one side is '-' => 3' overhang (0.5)
                           both sides have bases => internal (1.0)
  - Position penalty: +1 if reference 10th or 11th base differs
  - Length penalty: abs(len(query)-len(reference)) * 1.0
Total = 100 - 10 * (mismatch_penalty + position_penalty + length_penalty)
"""

import sys
import argparse

def read_fasta_seqs(file):
    """Return dict: name -> upper-case sequence."""
    seqs = {}
    cur_id = None
    cur_seq = []
    with open(file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if cur_id is not None:
                    seqs[cur_id] = ''.join(cur_seq).upper()
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line.upper())
        if cur_id is not None:
            seqs[cur_id] = ''.join(cur_seq).upper()
    return seqs

def revcomp(seq):
    """Reverse complement of a DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(comp.get(b, b) for b in reversed(seq))

def compute_penalties(qseq, sseq, qstart, qend, sstart, send, q_aln, s_aln):
    # ---------- Build global alignment (right‑align overhangs) ----------
    left_q = qstart - 1
    left_s = sstart - 1
    left_len = max(left_q, left_s)
    # Right‑align: pad left with gaps
    left_q_str = qseq[:left_q].rjust(left_len, '-') if left_q > 0 else '-' * left_len
    left_s_str = sseq[:left_s].rjust(left_len, '-') if left_s > 0 else '-' * left_len

    right_q = len(qseq) - qend
    right_s = len(sseq) - send
    right_len = max(right_q, right_s)
    right_q_str = qseq[qend:].ljust(right_len, '-') if right_q > 0 else '-' * right_len
    right_s_str = sseq[send:].ljust(right_len, '-') if right_s > 0 else '-' * right_len

    q_global = left_q_str + q_aln + right_q_str
    s_global = left_s_str + s_aln + right_s_str

    # ---------- Penalty computation (unchanged) ----------
    internal_pen = 0.0
    overhang5_pen = 0.0
    overhang3_pen = 0.0
    pos_pen = 0
    ref_pos = 0   # 1‑based on reference

    left_end = left_len
    internal_start = left_len
    internal_end = left_len + len(q_aln)
    right_start = internal_end

    # Left overhang
    for i in range(left_end):
        qb = q_global[i]
        sb = s_global[i]
        if sb != '-':
            ref_pos += 1
        if qb != sb:
            if qb == '-' or sb == '-':
                overhang5_pen += 1.0
            else:
                internal_pen += 1.0
        if (ref_pos == 10 or ref_pos == 11) and qb != sb:
            pos_pen += 1

    # Internal
    for i in range(internal_start, internal_end):
        qb = q_global[i]
        sb = s_global[i]
        if sb != '-':
            ref_pos += 1
        if qb != sb:
            internal_pen += 1.0
        if (ref_pos == 10 or ref_pos == 11) and qb != sb:
            pos_pen += 1

    # Right overhang
    for i in range(right_start, len(q_global)):
        qb = q_global[i]
        sb = s_global[i]
        if sb != '-':
            ref_pos += 1
        if qb != sb:
            if qb == '-' or sb == '-':
                overhang3_pen += 0.5
            else:
                internal_pen += 1.0
        if (ref_pos == 10 or ref_pos == 11) and qb != sb:
            pos_pen += 1

    mismatch_penalty = internal_pen + overhang5_pen + overhang3_pen
    length_penalty = abs(len(qseq) - len(sseq)) * 1.0
    return mismatch_penalty, pos_pen, length_penalty


def main():
    parser = argparse.ArgumentParser(description='Accurate miRNA alignment scoring')
    parser.add_argument('-i', '--input', required=True,
                        help='BLAST outfmt 6 with qseq sseq (13 columns)')
    parser.add_argument('-f', '--fasta', required=True,
                        help='FASTA with full sequences')
    parser.add_argument('-o', '--output', required=True,
                        help='Output score file')
    args = parser.parse_args()

    seqs = read_fasta_seqs(args.fasta)

    with open(args.input) as fin, open(args.output, 'w') as fout:
        fout.write("query\treference\tmismatch_penalty\tposition_penalty\tlength_penalty\ttotal_score\n")
        for line in fin:
            parts = line.strip().split('\t')
            if len(parts) < 13:
                continue
            qid, sid, _, _, _, _, qstart, qend, sstart, send, _, q_aln, s_aln = parts[:13]
            qstart, qend = int(qstart), int(qend)
            sstart, send = int(sstart), int(send)

            if qid not in seqs or sid not in seqs:
                sys.stderr.write(f"Warning: missing seq {qid} or {sid}\n")
                continue

            qseq = seqs[qid]
            sseq_orig = seqs[sid]

            # Handle reverse strand
            if sstart > send:
                sseq_full = revcomp(sseq_orig)
                slen = len(sseq_orig)
                new_sstart = slen - sstart + 1
                new_send = slen - send + 1
                sstart, send = new_sstart, new_send
            else:
                sseq_full = sseq_orig

            mismatch_pen, pos_pen, len_pen = compute_penalties(
                qseq, sseq_full, qstart, qend, sstart, send, q_aln, s_aln
            )

            total = 100.0 - 10.0 * (mismatch_pen + pos_pen + len_pen)
            fout.write(f"{qid}\t{sid}\t{mismatch_pen:.1f}\t{pos_pen}\t{len_pen:.1f}\t{total:.1f}\n")

if __name__ == '__main__':
    main()
