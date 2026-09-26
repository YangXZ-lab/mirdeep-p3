#!/usr/bin/env python3
"""
Script to unify miRNA names across multiple *basic-info files based on genomic
coordinates. It can operate in three modes:
1) With a reference file (-b/--basic): rename miRNAs in input files to match
   the reference names if their precursor and mature positions are within user-
   specified distances.
2) Without a reference file: merge all miRNAs from input files, cluster them by
   positional overlap, and assign continuous letter suffixes per family.
3) Common mode (-c/--common with -b): first use the reference to standardise
   names, then perform cross-file merging ONLY on entries that did NOT match
   the reference, preserving the reference-based names of matched entries.

Only families of the form MIR<digits> (e.g. MIR156, MIR319) are processed.
Families like MIRN1, MIRN2, or any non-numeric MIR* are left untouched.
"""

import argparse
import sys
import os
import re
from collections import defaultdict

# ------------------------------- data structures -------------------------------

class MiRNAEntry:
    """Holds information for one miRNA precursor line."""
    __slots__ = ('line', 'fields', 'family_str', 'family_num',
                 'chrom', 'strand', 'pre_start', 'pre_end',
                 'mat1_name', 'mat1_start', 'mat1_end',
                 'mat2_name', 'mat2_start', 'mat2_end',
                 'matched', 'is_standard')

    def __init__(self, line, fields):
        self.line = line
        self.fields = fields
        # family string e.g. 'MIR156', 'MIRN1'
        self.family_str = fields[2]
        # numeric part, e.g. '156' or 'N1' (used only for standard families)
        self.family_num = self.family_str[3:] if self.family_str.startswith('MIR') else self.family_str
        self.chrom = fields[4]
        self.strand = fields[5]
        self.pre_start = int(fields[6])
        self.pre_end = int(fields[7])
        # main mature
        self.mat1_name = fields[14]
        self.mat1_start = int(fields[15])
        self.mat1_end = int(fields[16])
        # star/passenger mature (may be absent)
        if len(fields) > 18 and fields[18].strip():
            self.mat2_name = fields[18]
            self.mat2_start = int(fields[19])
            self.mat2_end = int(fields[20])
        else:
            self.mat2_name = None
            self.mat2_start = None
            self.mat2_end = None
        self.matched = False  # whether it was matched to a reference entry
        # Only process families with a numeric suffix after 'MIR'
        self.is_standard = bool(re.match(r'^MIR\d+$', self.family_str))

    def get_prefix(self):
        """Return the current prefix from the first column (e.g. 'N710')."""
        return self.fields[0].split('-')[0] if '-' in self.fields[0] else self.fields[0]

    def set_names(self, new_full_prefix, new_family_num, letter):
        """Update all name fields (only called for standard families)."""
        self.fields[0] = f"{new_full_prefix}-MIR{new_family_num}{letter}"
        # main mature
        self.fields[14] = f"{new_full_prefix}-miR{new_family_num}{letter}"
        self.mat1_name = self.fields[14]
        # star mature if present
        if self.mat2_name is not None:
            self.fields[18] = f"{new_full_prefix}-miR{new_family_num}{letter}*"
            self.mat2_name = self.fields[18]

    def copy_names_from(self, other_entry):
        """Adopt all names from another entry (used when matched to reference)."""
        self.fields[0] = other_entry.fields[0]
        self.fields[14] = other_entry.fields[14]
        self.mat1_name = self.fields[14]
        if other_entry.mat2_name is not None and self.mat2_name is not None:
            self.fields[18] = other_entry.fields[18]
            self.mat2_name = self.fields[18]
        elif other_entry.mat2_name is not None and self.mat2_name is None:
            pass
        elif other_entry.mat2_name is None and self.mat2_name is not None:
            pass

    def get_main_mature_letter(self):
        """Extract the letter suffix from the main mature name (e.g. 'a' from 'miR156a')."""
        m = re.search(r'miR\d+([a-z]+)\*?$', self.mat1_name)
        return m.group(1) if m else ''

    def to_line(self):
        return '\t'.join(self.fields)


def parse_line(line):
    """Convert a raw line into a MiRNAEntry. Returns None if parsing fails."""
    line = line.strip()
    if not line:
        return None
    parts = line.split('\t')
    if len(parts) < 17:
        sys.stderr.write(f"Warning: line has fewer than 17 fields, skipped: {line[:80]}...\n")
        return None
    try:
        return MiRNAEntry(line, parts)
    except (IndexError, ValueError) as e:
        sys.stderr.write(f"Warning: skipping line due to parse error: {e}\n{line[:200]}")
        return None


# ------------------------------- matching logic -------------------------------

def mature_match(start1, end1, start2, end2, dist):
    return (abs(start1 - start2) <= dist and abs(end1 - end2) <= dist)

def entries_match(e1, e2, pre_dist, mat_dist):
    """Return True if two entries represent the same miRNA locus (both must be standard)."""
    # This function should only be called for standard entries; non-standard are never matched.
    if e1.family_str != e2.family_str:
        return False
    if e1.chrom != e2.chrom or e1.strand != e2.strand:
        return False
    if abs(e1.pre_start - e2.pre_start) > pre_dist:
        return False
    if abs(e1.pre_end - e2.pre_end) > pre_dist:
        return False
    if not mature_match(e1.mat1_start, e1.mat1_end,
                        e2.mat1_start, e2.mat1_end, mat_dist):
        return False
    return True


# ------------------------------- helper functions -------------------------------

def num_to_letters(n):
    """0->'a', 1->'b', ..., 25->'z', 26->'aa', ..."""
    result = ''
    while True:
        n, rem = divmod(n, 26)
        result = chr(ord('a') + rem) + result
        if n == 0:
            break
        n -= 1
    return result

def letter_to_num(letter):
    """Convert letter suffix back to 0-based index."""
    idx = 0
    for ch in letter:
        idx = idx * 26 + (ord(ch) - ord('a') + 1)
    return idx - 1

def max_letter_index_from_entries(entries, family_str):
    """Return the maximum 0-based letter index for a given family among a list of entries."""
    max_idx = -1
    for e in entries:
        if e.family_str != family_str:
            continue
        letter = e.get_main_mature_letter()
        if not letter:
            continue
        idx = letter_to_num(letter)
        if idx > max_idx:
            max_idx = idx
    return max_idx

def chrom_sort_key(chrom):
    """Natural sort key for chromosome names (e.g. Chr03 < Chr04 < Chr11)."""
    m = re.match(r'^([A-Za-z]+)(\d+)$', chrom)
    if m:
        return (m.group(1), int(m.group(2)))
    return (chrom, 0)


# ------------------------------- Union-Find for clustering -------------------------------

class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, x, y):
        rx = self.find(x)
        ry = self.find(y)
        if rx == ry:
            return
        if self.rank[rx] < self.rank[ry]:
            self.parent[rx] = ry
        elif self.rank[rx] > self.rank[ry]:
            self.parent[ry] = rx
        else:
            self.parent[ry] = rx
            self.rank[rx] += 1


# ------------------------------- cross-file merging of unmatched entries -------------------------------

def merge_unmatched_entries(unmatched_list, ref_entries, pre_dist, mat_dist):
    """
    Given a flat list of standard MiRNAEntry objects that were NOT matched to the
    reference, cluster them by positional overlap (within the same family) and
    rename them with continuous letters starting after the reference's maximum
    letter for that family.
    """
    if not unmatched_list:
        return

    by_family = defaultdict(list)
    for e in unmatched_list:
        by_family[e.family_str].append(e)

    for fam, entries in by_family.items():
        n = len(entries)
        uf = UnionFind(n)
        for i in range(n):
            for j in range(i + 1, n):
                if entries_match(entries[i], entries[j], pre_dist, mat_dist):
                    uf.union(i, j)

        clusters = defaultdict(list)
        for idx, e in enumerate(entries):
            root = uf.find(idx)
            clusters[root].append(e)

        rep_list = []
        for root, members in clusters.items():
            rep = min(members, key=lambda x: (chrom_sort_key(x.chrom), x.pre_start, x.mat1_start))
            rep_list.append((rep, members))
        rep_list.sort(key=lambda x: (chrom_sort_key(x[0].chrom), x[0].pre_start, x[0].mat1_start))

        # Starting letter index: reference max + 1
        start_idx = max_letter_index_from_entries(ref_entries, fam) + 1
        # Determine prefix: try to get from reference, otherwise from first entry
        ref_prefix = None
        for ref_e in ref_entries:
            if ref_e.family_str == fam:
                ref_prefix = ref_e.get_prefix()
                break
        default_prefix = entries[0].get_prefix() if entries else ''

        for i, (_, members) in enumerate(rep_list):
            letter = num_to_letters(start_idx + i)
            prefix = ref_prefix if ref_prefix is not None else default_prefix
            for e in members:
                e.set_names(prefix, e.family_num, letter)


# ------------------------------- mode 1: with reference -------------------------------

def process_with_reference(input_files, ref_file, pre_dist, mat_dist, output_files):
    # Load reference, but keep only standard families
    ref_entries = []
    with open(ref_file, 'r') as fh:
        for line in fh:
            e = parse_line(line)
            if e and e.is_standard:
                ref_entries.append(e)
    if not ref_entries:
        sys.stderr.write("Warning: reference file contains no standard miRNA entries.\n")

    for in_path, out_path in zip(input_files, output_files):
        in_entries = []
        with open(in_path, 'r') as fh:
            for line in fh:
                e = parse_line(line)
                if e:
                    in_entries.append(e)
                else:
                    in_entries.append(line)

        unmatched = []
        for entry in in_entries:
            if not isinstance(entry, MiRNAEntry) or not entry.is_standard:
                continue
            matched_ref = None
            for ref_e in ref_entries:
                if entries_match(entry, ref_e, pre_dist, mat_dist):
                    matched_ref = ref_e
                    break
            if matched_ref:
                entry.copy_names_from(matched_ref)
                entry.matched = True
            else:
                entry.matched = False
                unmatched.append(entry)

        # Temporary naming for unmatched standard entries (will not be used if -c is active)
        by_family = defaultdict(list)
        for e in unmatched:
            by_family[e.family_str].append(e)
        for fam, fam_entries in by_family.items():
            max_idx = max_letter_index_from_entries(ref_entries, fam)
            fam_entries.sort(key=lambda x: (chrom_sort_key(x.chrom), x.pre_start))
            ref_prefix = None
            for ref_e in ref_entries:
                if ref_e.family_str == fam:
                    ref_prefix = ref_e.get_prefix()
                    break
            for i, e in enumerate(fam_entries):
                new_letter = num_to_letters(max_idx + 1 + i)
                prefix = ref_prefix if ref_prefix is not None else e.get_prefix()
                e.set_names(prefix, e.family_num, new_letter)

        os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
        with open(out_path, 'w') as out_fh:
            for item in in_entries:
                if isinstance(item, MiRNAEntry):
                    out_fh.write(item.to_line() + '\n')
                else:
                    out_fh.write(item if isinstance(item, str) else str(item))


# ------------------------------- mode 2: merge without reference -------------------------------

def process_merge_only(input_files, pre_dist, mat_dist, output_files):
    file_entries = []
    all_entries = []  # only standard entries

    for in_path in input_files:
        entries = []
        with open(in_path, 'r') as fh:
            for line in fh:
                e = parse_line(line)
                if e:
                    entries.append(e)
                    if e.is_standard:
                        all_entries.append(e)
                else:
                    entries.append(line)
        file_entries.append(entries)

    if all_entries:
        n = len(all_entries)
        uf = UnionFind(n)
        for i in range(n):
            for j in range(i + 1, n):
                if entries_match(all_entries[i], all_entries[j], pre_dist, mat_dist):
                    uf.union(i, j)

        clusters = defaultdict(list)
        for idx, e in enumerate(all_entries):
            root = uf.find(idx)
            clusters[root].append(e)

        family_clusters = defaultdict(list)
        for root, members in clusters.items():
            rep = min(members, key=lambda x: (chrom_sort_key(x.chrom), x.pre_start, x.mat1_start))
            family_clusters[rep.family_str].append((rep, members))

        for fam, cluster_data in family_clusters.items():
            cluster_data.sort(key=lambda x: (chrom_sort_key(x[0].chrom),
                                             x[0].pre_start,
                                             x[0].mat1_start))
            for idx, (rep, members) in enumerate(cluster_data):
                letter = num_to_letters(idx)
                for e in members:
                    prefix = e.get_prefix()
                    e.set_names(prefix, e.family_num, letter)

    for entries, out_path in zip(file_entries, output_files):
        os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
        with open(out_path, 'w') as out_fh:
            for item in entries:
                if isinstance(item, MiRNAEntry):
                    out_fh.write(item.to_line() + '\n')
                else:
                    out_fh.write(item if isinstance(item, str) else str(item))


# ------------------------------- mode 3: common (reference + merge unmatched) -------------------------------

def process_common(input_files, ref_file, pre_dist, mat_dist, outputb_files, output_files):
    # Load reference, standard only
    ref_entries = []
    with open(ref_file, 'r') as fh:
        for line in fh:
            e = parse_line(line)
            if e and e.is_standard:
                ref_entries.append(e)

    intermediate_file_entries = []
    all_unmatched = []   # will only contain standard unmatched entries

    # ---- Step 1: reference-based renaming ----
    for in_path in input_files:
        in_entries = []
        unmatched_per_file = []
        with open(in_path, 'r') as fh:
            for line in fh:
                e = parse_line(line)
                if e:
                    in_entries.append(e)
                else:
                    in_entries.append(line)

        for entry in in_entries:
            if not isinstance(entry, MiRNAEntry) or not entry.is_standard:
                continue
            matched_ref = None
            for ref_e in ref_entries:
                if entries_match(entry, ref_e, pre_dist, mat_dist):
                    matched_ref = ref_e
                    break
            if matched_ref:
                entry.copy_names_from(matched_ref)
                entry.matched = True
            else:
                entry.matched = False
                unmatched_per_file.append(entry)

        # Temporary per-file naming for unmatched standard entries (will be overridden globally)
        by_family = defaultdict(list)
        for e in unmatched_per_file:
            by_family[e.family_str].append(e)
        for fam, fam_entries in by_family.items():
            max_idx = max_letter_index_from_entries(ref_entries, fam)
            fam_entries.sort(key=lambda x: (chrom_sort_key(x.chrom), x.pre_start))
            ref_prefix = None
            for ref_e in ref_entries:
                if ref_e.family_str == fam:
                    ref_prefix = ref_e.get_prefix()
                    break
            for i, e in enumerate(fam_entries):
                new_letter = num_to_letters(max_idx + 1 + i)
                prefix = ref_prefix if ref_prefix is not None else e.get_prefix()
                e.set_names(prefix, e.family_num, new_letter)

        intermediate_file_entries.append(in_entries)
        all_unmatched.extend(unmatched_per_file)

    # ---- Optional output of intermediate files ----
    if outputb_files:
        for in_entries, out_path in zip(intermediate_file_entries, outputb_files):
            os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
            with open(out_path, 'w') as out_fh:
                for item in in_entries:
                    if isinstance(item, MiRNAEntry):
                        out_fh.write(item.to_line() + '\n')
                    else:
                        out_fh.write(item if isinstance(item, str) else str(item))

    # ---- Step 2: cross-file merging of unmatched standard entries only ----
    if all_unmatched:
        merge_unmatched_entries(all_unmatched, ref_entries, pre_dist, mat_dist)

    # ---- Final output ----
    if output_files:
        for in_entries, out_path in zip(intermediate_file_entries, output_files):
            os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
            with open(out_path, 'w') as out_fh:
                for item in in_entries:
                    if isinstance(item, MiRNAEntry):
                        out_fh.write(item.to_line() + '\n')
                    else:
                        out_fh.write(item if isinstance(item, str) else str(item))


# ------------------------------- command line interface -------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Unify miRNA names across multiple basic-info files based on genomic positions.')
    parser.add_argument('-i', '--input', required=True,
                        help='Comma-separated list of input files')
    parser.add_argument('-b', '--basic', default=None,
                        help='Reference basic-info file (e.g. Ath-basic-info)')
    parser.add_argument('--mature-distance', type=int, default=1,
                        help='Allowed mismatch distance for mature start/end (default: 1)')
    parser.add_argument('--pre-distance', type=int, default=2,
                        help='Allowed mismatch distance for precursor start/end (default: 2)')
    parser.add_argument('-ob', '--outputb', default=None,
                        help='Comma-separated output files after reference-based renaming')
    parser.add_argument('-o', '--output', default=None,
                        help='Comma-separated output files for final (or only) renaming')
    parser.add_argument('-c', '--common', action='store_true',
                        help='Combine reference-based renaming with cross-file merging of unmatched entries (requires -b)')

    args = parser.parse_args()

    input_files = [f.strip() for f in args.input.split(',') if f.strip()]
    if not input_files:
        sys.exit("Error: at least one input file must be provided via -i/--input.")

    for f in input_files:
        if not os.path.isfile(f):
            sys.exit(f"Error: input file '{f}' does not exist.")
    if args.basic and not os.path.isfile(args.basic):
        sys.exit(f"Error: reference file '{args.basic}' does not exist.")

    if args.common and not args.basic:
        sys.exit("Error: --common/-c requires -b/--basic to be specified.")

    # Determine mode and validate outputs
    if args.basic and not args.common:
        # Mode 1: reference only
        if not args.outputb:
            sys.exit("Error: -ob/--outputb must be specified when using -b without -c.")
        outputb_list = [f.strip() for f in args.outputb.split(',') if f.strip()]
        if len(outputb_list) != len(input_files):
            sys.exit(f"Error: number of outputb files ({len(outputb_list)}) must match input files ({len(input_files)}).")
        process_with_reference(input_files, args.basic,
                               args.pre_distance, args.mature_distance,
                               outputb_list)
    elif args.common:
        # Mode 3: common
        if not args.outputb and not args.output:
            sys.exit("Error: at least one of -ob/--outputb or -o/--output must be specified in common mode.")
        outputb_list = []
        output_list = []
        if args.outputb:
            outputb_list = [f.strip() for f in args.outputb.split(',') if f.strip()]
            if len(outputb_list) != len(input_files):
                sys.exit(f"Error: number of outputb files ({len(outputb_list)}) must match input files ({len(input_files)}).")
        if args.output:
            output_list = [f.strip() for f in args.output.split(',') if f.strip()]
            if len(output_list) != len(input_files):
                sys.exit(f"Error: number of output files ({len(output_list)}) must match input files ({len(input_files)}).")
        process_common(input_files, args.basic,
                       args.pre_distance, args.mature_distance,
                       outputb_list if outputb_list else None,
                       output_list if output_list else None)
    else:
        # Mode 2: merge only (no basic)
        if not args.output:
            sys.exit("Error: -o/--output must be specified when -b is absent.")
        output_list = [f.strip() for f in args.output.split(',') if f.strip()]
        if len(output_list) != len(input_files):
            sys.exit(f"Error: number of output files ({len(output_list)}) must match input files ({len(input_files)}).")
        process_merge_only(input_files,
                           args.pre_distance, args.mature_distance,
                           output_list)


if __name__ == '__main__':
    main()