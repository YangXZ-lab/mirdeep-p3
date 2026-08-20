#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.


"""
filter_alignments.py - Filter BLAST-parsed alignment results

This script parses a BLAST-parsed alignment file and filters alignments based on criteria.
By default, it prints lines where queries align perfectly to subjects.
Options:
  -a: Allow mismatched nucleotides in the 3' end of the query
  -b: Print filtered reads in fasta format
  -c: Filter sequences mapping more than specified number of times
"""

import sys
import os
import re
from collections import defaultdict
from typing import Dict, List, Tuple, Set
import argparse

class Alignment:
    """Represents a single alignment record"""
    
    def __init__(self, line: str):
        """Parse a BLAST-parsed line"""
        self.line = line.strip()
        parts = self.line.split('\t')
        
        # BLAST format has 10 fields, not 12
        if len(parts) < 10:
            raise ValueError(f"Invalid BLAST line: expected at least 10 fields, got {len(parts)}")
        
        self.query_id = parts[0]
        self.query_length = int(parts[1])
        self.query_range = parts[2]
        self.subject_id = parts[3]
        self.subject_length = int(parts[4])
        self.subject_range = parts[5]
        self.e_value = parts[6]
        self.identity = float(parts[7])
        self.bitscore = float(parts[8])
        # The last field may contain strand info
        self.strand_info = parts[9] if len(parts) > 9 else ""
        
        # Parse query range
        self.query_beg, self.query_end = self._parse_range(self.query_range)
        
        # Parse subject range
        self.subject_beg, self.subject_end = self._parse_range(self.subject_range)
        
        # Calculate truncation
        self.trunc_beg = self.query_beg - 1  # uncovered at 5' end
        self.trunc_end = self.query_length - self.query_end  # uncovered at 3' end
        
        # Determine if it's a perfect alignment based on strand
        self.is_perfect = (self.identity == 1.0)
    
    def _parse_range(self, range_str: str) -> Tuple[int, int]:
        """Parse range string like '1..21' into (start, end)"""
        match = re.match(r'(\d+)\.\.(\d+)', range_str)
        if not match:
            raise ValueError(f"Invalid range format: {range_str}")
        return int(match.group(1)), int(match.group(2))
    
    def __repr__(self):
        return f"Alignment(query={self.query_id}, subject={self.subject_id}, identity={self.identity}, trunc_end={self.trunc_end})"
    
    def __str__(self):
        return self.line

class AlignmentFilter:
    """Main class for filtering alignments"""
    
    def __init__(self, trunc_end_max: int = 0, max_mappings: int = None):
        """
        Initialize filter
        
        Args:
            trunc_end_max: Maximum allowed uncovered nucleotides at 3' end
            max_mappings: Maximum allowed mappings per query (None for no limit)
        """
        self.trunc_end_max = trunc_end_max
        self.max_mappings = max_mappings
        
        # Data structures
        self.alignments_by_query = defaultdict(list)  # query -> list of alignments
        self.best_trunc_by_query = {}  # query -> best trunc_end value
        self.query_sequences = {}  # query -> sequence (if -b option used)
    
    def parse_blast_file(self, blast_file: str):
        """Parse BLAST-parsed file and store alignments"""
        try:
            total_lines = 0
            valid_lines = 0
            with open(blast_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:
                        continue
                    
                    total_lines += 1
                    
                    try:
                        alignment = Alignment(line)
                        
                        # Check identity (must be 1.0 for perfect alignment)
                        if not alignment.is_perfect:
                            continue
                        
                        # Check 5' truncation (must be 0)
                        if alignment.trunc_beg != 0:
                            continue
                        
                        # Check 3' truncation (must be <= trunc_end_max)
                        if alignment.trunc_end > self.trunc_end_max:
                            continue
                        
                        # Store alignment
                        query = alignment.query_id
                        self.alignments_by_query[query].append(alignment)
                        
                        # Update best truncation for this query
                        if (query not in self.best_trunc_by_query or 
                            alignment.trunc_end < self.best_trunc_by_query[query]):
                            self.best_trunc_by_query[query] = alignment.trunc_end
                        
                        valid_lines += 1
                            
                    except ValueError as e:
                        sys.stderr.write(f"Warning: Skipping line {line_num}: {e}\n")
                        continue
            
            print(f"Total lines read: {total_lines:,}")
            print(f"Valid alignments parsed: {valid_lines:,}")
            print(f"Found {len(self.alignments_by_query):,} unique queries")
            
            # Debug: show first few queries and their alignment counts
            if len(self.alignments_by_query) > 0:
                print(f"Sample query counts (first 5):")
                for i, (query, aligns) in enumerate(list(self.alignments_by_query.items())[:5]):
                    print(f"  {query}: {len(aligns)} alignments")
            
        except Exception as e:
            sys.stderr.write(f"Error reading BLAST file {blast_file}: {e}\n")
            sys.exit(1)
    
    def parse_fasta_file(self, fasta_file: str):
        """Parse FASTA file to get query sequences"""
        try:
            with open(fasta_file, 'r') as f:
                current_id = None
                current_seq = []
                
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    if line.startswith('>'):
                        # Save previous sequence
                        if current_id is not None:
                            self.query_sequences[current_id] = ''.join(current_seq)
                        
                        # Start new sequence
                        current_id = line[1:].split()[0]  # Get ID (first word after >)
                        current_seq = []
                    else:
                        current_seq.append(line)
                
                # Save last sequence
                if current_id is not None:
                    self.query_sequences[current_id] = ''.join(current_seq)
            
            print(f"Parsed {len(self.query_sequences):,} sequences from FASTA file")
            
        except Exception as e:
            sys.stderr.write(f"Error reading FASTA file {fasta_file}: {e}\n")
            sys.exit(1)
    
    def filter_by_mapping_count(self, query: str) -> bool:
        """
        Check if a query passes the -c filter (max mappings)
        
        Returns:
            True if query passes filter, False otherwise
        """
        if self.max_mappings is None:
            return True
        
        # Get best truncation for this query
        best_trunc = self.best_trunc_by_query.get(query)
        if best_trunc is None:
            return False
        
        # Count alignments with best truncation
        count = 0
        for alignment in self.alignments_by_query[query]:
            if alignment.trunc_end == best_trunc:
                count += 1
        
        return count <= self.max_mappings
    
    def output_blast_format(self, output_file: str = None):
        """Output filtered alignments in BLAST format"""
        output_handle = sys.stdout
        if output_file:
            try:
                output_handle = open(output_file, 'w')
                print(f"Outputting filtered alignments to: {output_file}")
            except Exception as e:
                sys.stderr.write(f"Error opening output file {output_file}: {e}\n")
                sys.exit(1)
        else:
            print("Outputting filtered alignments to stdout...")
        
        # Sort queries for consistent output
        sorted_queries = sorted(self.alignments_by_query.keys())
        
        output_count = 0
        for query in sorted_queries:
            # Check mapping count filter
            if not self.filter_by_mapping_count(query):
                continue
            
            # Get best truncation for this query
            best_trunc = self.best_trunc_by_query[query]
            
            # Output all alignments with best truncation
            for alignment in self.alignments_by_query[query]:
                if alignment.trunc_end == best_trunc:
                    output_handle.write(str(alignment) + "\n")
                    output_count += 1
        
        if output_file:
            output_handle.close()
        
        print(f"Output {output_count:,} alignments")
        return output_count
    
    def output_fasta_format(self, output_file: str = None):
        """Output filtered queries in FASTA format"""
        output_handle = sys.stdout
        if output_file:
            try:
                output_handle = open(output_file, 'w')
                print(f"Outputting filtered sequences to: {output_file}")
            except Exception as e:
                sys.stderr.write(f"Error opening output file {output_file}: {e}\n")
                sys.exit(1)
        else:
            print("Outputting filtered sequences to stdout...")
        
        # Sort queries for consistent output
        sorted_queries = sorted(self.alignments_by_query.keys())
        
        output_count = 0
        for query in sorted_queries:
            # Check mapping count filter
            if not self.filter_by_mapping_count(query):
                continue
            
            # Get sequence
            sequence = self.query_sequences.get(query)
            if not sequence:
                sys.stderr.write(f"Warning: Sequence not found for query {query}\n")
                continue
            
            # Get best truncation for this query
            best_trunc = self.best_trunc_by_query[query]
            
            # Truncate sequence if needed
            if best_trunc > 0:
                truncated_seq = sequence[:-best_trunc]  # Remove from 3' end
                suffix = f"_t{best_trunc}"
            else:
                truncated_seq = sequence
                suffix = ""
            
            # Output in FASTA format
            output_handle.write(f">{query}{suffix}\n")
            output_handle.write(f"{truncated_seq}\n")
            output_count += 1
        
        if output_file:
            output_handle.close()
        
        print(f"Output {output_count:,} sequences")
        return output_count
    
    def get_filtering_stats(self) -> Dict:
        """Get statistics about filtering"""
        stats = {
            'total_queries': len(self.alignments_by_query),
            'total_alignments': sum(len(v) for v in self.alignments_by_query.values()),
        }
        
        # Count queries that pass -c filter
        if self.max_mappings is not None:
            passing_queries = sum(1 for q in self.alignments_by_query.keys() 
                                 if self.filter_by_mapping_count(q))
            stats['passing_queries'] = passing_queries
            stats['filtered_queries'] = stats['total_queries'] - passing_queries
        
        return stats

def main():
    parser = argparse.ArgumentParser(
        description='Filter BLAST-parsed alignment results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Input format (BLAST-parsed, tab-separated):
  query_id, query_length, query_range, subject_id, subject_length,
  subject_range, e-value, identity, bitscore, strand_info
  
Filtering criteria:
  1. Identity must be 1.0 (perfect alignment)
  2. No uncovered nucleotides at 5' end (trunc_beg = 0)
  3. Limited uncovered nucleotides at 3' end (controlled by -a)
  4. Maximum mapping count (controlled by -c)
  
Examples:
  # Default: output perfect alignments
  python filter_alignments.py input.bst
  
  # Allow up to 2 uncovered nucleotides at 3' end
  python filter_alignments.py input.bst -a 2
  
  # Filter queries mapping more than 10 times
  python filter_alignments.py input.bst -c 10 -o output.bst
  
  # Output filtered sequences in FASTA format
  python filter_alignments.py input.bst -b sequences.fa -o output.fa
        """
    )
    
    parser.add_argument('blast_file',
                       help='BLAST-parsed alignment file')
    
    parser.add_argument('-a', '--allow-trunc-end', type=int, default=0,
                       help='Maximum uncovered nucleotides allowed at 3\' end (default: 0)')
    
    parser.add_argument('-b', '--fasta-file',
                       help='Input FASTA file (output filtered sequences instead of alignments)')
    
    parser.add_argument('-c', '--max-mappings', type=int,
                       help='Maximum number of mappings allowed per query')
    
    parser.add_argument('-o', '--output',
                       help='Output file (default: stdout)')
    
    parser.add_argument('--stats', action='store_true',
                       help='Show detailed filtering statistics')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("FILTERING BLAST ALIGNMENTS")
    print("=" * 60)
    print(f"Input file: {args.blast_file}")
    print(f"Parameters:")
    print(f"  Allow 3' truncation: {args.allow_trunc_end} nt")
    if args.max_mappings:
        print(f"  Max mappings per query: {args.max_mappings}")
    if args.fasta_file:
        print(f"  FASTA file for output: {args.fasta_file}")
    if args.output:
        print(f"  Output file: {args.output}")
    else:
        print(f"  Output: stdout")
    print("-" * 60)
    
    # Initialize filter
    filter_obj = AlignmentFilter(
        trunc_end_max=args.allow_trunc_end,
        max_mappings=args.max_mappings
    )
    
    # Parse BLAST file
    print("1. Parsing BLAST file...")
    filter_obj.parse_blast_file(args.blast_file)
    
    # Parse FASTA file if -b option used
    if args.fasta_file:
        print("2. Parsing FASTA file...")
        filter_obj.parse_fasta_file(args.fasta_file)
    
    # Output results
    print("3. Filtering and outputting results...")
    if args.fasta_file:
        output_count = filter_obj.output_fasta_format(args.output)
    else:
        output_count = filter_obj.output_blast_format(args.output)
    
    # Show statistics
    if args.stats:
        print("\n" + "=" * 60)
        print("FILTERING STATISTICS")
        print("=" * 60)
        
        stats = filter_obj.get_filtering_stats()
        
        print(f"Total queries: {stats['total_queries']:,}")
        print(f"Total alignments: {stats['total_alignments']:,}")
        
        if 'passing_queries' in stats:
            print(f"\nWith -c {args.max_mappings} filter:")
            print(f"  Queries passing filter: {stats['passing_queries']:,}")
            print(f"  Queries filtered out: {stats['filtered_queries']:,}")
            print(f"  Filter rate: {(stats['filtered_queries']/stats['total_queries']*100):.1f}%")
        
        print(f"\nFinal output count: {output_count:,}")
        
        # Truncation distribution
        trunc_counts = defaultdict(int)
        for query, best_trunc in filter_obj.best_trunc_by_query.items():
            trunc_counts[best_trunc] += 1
        
        print(f"\nTruncation distribution (3' end):")
        for trunc_len in sorted(trunc_counts.keys()):
            count = trunc_counts[trunc_len]
            percentage = (count / stats['total_queries']) * 100
            print(f"  {trunc_len} nt uncovered: {count:,} queries ({percentage:.1f}%)")
    
    print("\n" + "=" * 60)
    print("FILTERING COMPLETE")
    print("=" * 60)

if __name__ == "__main__":
    main()
