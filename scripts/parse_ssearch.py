#!/usr/bin/env python

"""parse_ssearch.py. From the default output of ssearch36, this script parses each alignment to a .tsv format that can then be used in the companion 
script parse_mirna_targets.py
Developed and tested with python 3.9.0

NOTES:
It is not entirely clear which ssearch36 parameters are used by psRNATarget, but using the following ssearch command followed by the execution of parse_ssearch.py
and then followed by the execution of this script with default paramemeters managed to reproduced the output of psRNATarget.

ssearch36 -f -8 -g -3 -E 10000 -T 8 -b 200 -r +4/-3 -n -U -W 10 -N 20000 -i <input_mirna.fasta> <reference_genome.fasta> > <output_file>
If you provide an already rev-complemented fasta file, you can omit the -i argument.

Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
Modified by: Jiawen zhao
"""

import argparse
import os
import sys
import re
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def parse_command_line_arguments():
    parser = argparse.ArgumentParser(description='Convert SSEARCH36 default output to a .tsv format output')
    parser.add_argument('-i', '--infile', required=False, help='Input file (i.e. output of ssearch36). This argument is optional as the output of ssearch36 can be piped directly to this script as well.', type=argparse.FileType('r'))
    parser.add_argument('--rev', default=False, action=argparse.BooleanOptionalAction, help='Reverse')
    parser.add_argument('--verbose', default=False, action=argparse.BooleanOptionalAction, help='Verbose output')
    args = parser.parse_args()
    return args

def main(arguments):
    args = parse_command_line_arguments()
    verbose = args.verbose
    rev = args.rev
    
    # Dictionary to track already output alignments (optional warning instead of error)
    seen = {}
    curr_query = ""
    curr_target = ""
    curr_aln = ""          # not used explicitly but kept for clarity
    curr_query_length = 0
    counter_match_nucl_string = 0
    curr_start = 0
    curr_end = 0
    start = 0
    end = 0
    hsp = 0
    q_end = 0
    q_start = 0
    query_str = ""
    aln_str = ""
    subject_str = ""
    i = 0
    strand = ""
    
    if args.infile:
        infile = os.path.abspath(args.infile.name)
        fhand = open(infile, 'r')
    else:
        fhand = sys.stdin

    for line in fhand:
        line = line.rstrip()

        if line.startswith("#"):
            continue

        # ---- Relaxed query line matching: only capture the query name ----
        match = re.match(r"\s*\d+>>>(\S+)", line)
        if match:
            if verbose:
                print("\n---------------------------------------------", file=sys.stderr)
                print("curr_query: " + match.group(1), file=sys.stderr)
            # Reset state for a new query
            curr_query = match.group(1)
            counter_match_nucl_string = 0
            curr_start = 0
            curr_end = 0
            curr_query_length = 0  # will be filled from alignment later
            i = 0
            continue

        # ---- Relaxed target line matching: only capture the first word ----
        match = re.match(r"^>>(\S+)", line)
        if match:
            i += 1
            curr_target = match.group(1) + "_" + str(i)
            continue

        # ---- Smith-Waterman score line (unchanged, but relies on curr_query_length) ----
        match = re.match(r"^Smith-Waterman.*\((\d+)-(\d+):(\d+)-(\d+)\)$", line)
        if match:
            # At this point curr_query_length must be set (from query_str)
            diff_start = int(match.group(1)) - 1
            diff_end = int(curr_query_length)
            start = int(match.group(3)) - diff_start
            end = int(match.group(4)) + diff_end
            hsp = int(match.group(4)) - int(match.group(3))
            if rev is True:
                q_start = int(match.group(1))
                q_end = int(match.group(2))
                strand = "-"
            else:
                # If fwd, qstart is in reverse orientation.
                q_start = int(match.group(2))
                q_end = int(match.group(1))
                strand = "+"
            continue

        # ---- Ignore separator line between hits of the same contig ----
        match = re.match(r"^>--$", line)
        if match:
            i += 1
            continue

        # ---- Query nucleotide string (first line of alignment block) ----
        match = re.match(r"^\S+\s+([ACGTU-]*)\s*$", line)
        if match and counter_match_nucl_string == 0:
            # Found a query nucl alignment string.
            curr_start = match.span(1)[0]
            curr_end = match.span(1)[1]
            curr_str = line[curr_start:curr_end]
            query_str = curr_str

            # *** NEW: derive query length from the aligned sequence (remove gaps) ***
            curr_query_length = len(query_str.replace('-', '').replace(' ', ''))

            # Record that we have seen this query-target pair (optional warning)
            pair_key = curr_target + "_" + str(i)
            if curr_query not in seen:
                seen[curr_query] = set()
            if pair_key in seen[curr_query]:
                if verbose:
                    print(f"Warning: duplicate alignment for {curr_query} {pair_key}", file=sys.stderr)
            else:
                seen[curr_query].add(pair_key)

            counter_match_nucl_string = 1
            continue

        # ---- Match line (dots/colons indicating identity/similarity) ----
        match = re.match(r"^\s+[\.\:]", line)
        if match:
            if counter_match_nucl_string == 1:
                # Found a match string.
                # Extract substring at previously found positions.
                curr_str = line[curr_start:curr_end]
                aln_str = curr_str
                continue

        # ---- Subject nucleotide string (last line of alignment block) ----
        match = re.match(r"^\S+\s+[ACGTBDHUYNVWRMKS-]*\s*$", line)
        if match and counter_match_nucl_string == 1:
            # Extract substring at previously found positions.
            curr_str = line[curr_start:curr_end]
            subject_str = curr_str

            # Safety check: make sure we have a valid query and target
            if not curr_query or not curr_target:
                if verbose:
                    print(f"Warning: incomplete alignment block, skipping", file=sys.stderr)
                counter_match_nucl_string = 0
                continue

            # Optional duplicate check (no longer raises an exception)
            pair_key = curr_target + "_" + str(i)
            if curr_query in seen and pair_key not in seen[curr_query]:
                if verbose:
                    print(f"Warning: unexpected alignment order for {curr_query} {pair_key}", file=sys.stderr)

            # Output the parsed alignment
            print(str(curr_query) + "\t" + str(curr_target) + "\t" + str(aln_str) + "\t" + str(query_str) + "\t" + str(subject_str) + "\t" + str(q_start) + "\t" + str(q_end) + "\t" + str(start) + "\t" + str(end) + "\t" + str(strand) + "\t" + str(hsp), file=sys.stdout)
            
            counter_match_nucl_string = 0
            continue

    fhand.close()

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
