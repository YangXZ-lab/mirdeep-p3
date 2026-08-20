#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
excise_candidate.py - Excise potential miRNA precursor sequences from genome

This script excises potential miRNA precursor sequences from a genome
using the positions of aligned reads as guidelines. It takes a genome
FASTA file, a BLAST-parsed alignment file, and a precursor length limit,
then outputs candidate precursor sequences in FASTA format.
"""

import sys
import os
import re
from collections import defaultdict
from typing import Dict, List, Tuple, Set
import argparse

def parse_fasta(fasta_file: str) -> Dict[str, str]:
    """
    Parse FASTA file and return dictionary of sequences.
    
    Args:
        fasta_file: Path to FASTA file
        
    Returns:
        Dictionary mapping sequence ID to sequence
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('>'):
                    # Save previous sequence
                    if current_id is not None:
                        sequences[current_id] = ''.join(current_seq)
                    
                    # Start new sequence
                    current_id = line[1:].split()[0]  # Get ID (first word after >)
                    current_seq = []
                else:
                    current_seq.append(line.upper())
            
            # Save last sequence
            if current_id is not None:
                sequences[current_id] = ''.join(current_seq)
        
        print(f"Parsed {len(sequences)} sequences from {fasta_file}")
        return sequences
        
    except Exception as e:
        sys.stderr.write(f"Error: Failed to parse FASTA file {fasta_file}: {e}\n")
        sys.exit(1)

def parse_blast_file(blast_file: str) -> List[Dict]:
    """
    Parse BLAST-parsed alignment file.
    
    Args:
        blast_file: Path to BLAST-parsed file
        
    Returns:
        List of alignment dictionaries
    """
    alignments = []
    
    try:
        with open(blast_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) < 10:
                    sys.stderr.write(f"Warning: Line {line_num} has insufficient fields: {line}\n")
                    continue
                
                try:
                    query_id = parts[0]
                    query_length = int(parts[1])
                    query_range = parts[2]
                    subject_id = parts[3]
                    subject_length = int(parts[4])
                    subject_range = parts[5]
                    e_value = parts[6]
                    identity = float(parts[7])
                    bitscore = float(parts[8])
                    strand_info = parts[9]
                    
                    # Parse ranges
                    query_match = re.match(r'(\d+)\.\.(\d+)', query_range)
                    subject_match = re.match(r'(\d+)\.\.(\d+)', subject_range)
                    
                    if not query_match or not subject_match:
                        sys.stderr.write(f"Warning: Invalid range format on line {line_num}\n")
                        continue
                    
                    query_beg = int(query_match.group(1))
                    query_end = int(query_match.group(2))
                    subject_beg = int(subject_match.group(1))
                    subject_end = int(subject_match.group(2))
                    
                    # Determine strand from strand_info
                    strand = '+'
                    if 'Minus' in strand_info or '-' in strand_info:
                        strand = '-'
                    
                    alignments.append({
                        'query': query_id,
                        'query_beg': query_beg,
                        'query_end': query_end,
                        'query_length': query_length,
                        'subject': subject_id,
                        'subject_beg': subject_beg,
                        'subject_end': subject_end,
                        'subject_length': subject_length,
                        'strand': strand,
                        'line': line
                    })
                    
                except (ValueError, IndexError) as e:
                    sys.stderr.write(f"Warning: Error parsing line {line_num}: {e}\n")
                    continue
        
        print(f"Parsed {len(alignments)} alignments from {blast_file}")
        return alignments
        
    except Exception as e:
        sys.stderr.write(f"Error: Failed to read BLAST file {blast_file}: {e}\n")
        sys.exit(1)

def revcom(sequence: str) -> str:
    """
    Return reverse complement of DNA sequence.
    
    Args:
        sequence: DNA sequence
        
    Returns:
        Reverse complement sequence
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 
                  'N': 'N', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))

def merge_features(features: Dict, alignment: Dict):
    """
    Merge alignment with existing features if they overlap or are close.
    
    Args:
        features: Dictionary containing current features
        alignment: New alignment to merge
    """
    subject = alignment['subject']
    strand = alignment['strand']
    subject_beg = alignment['subject_beg']
    subject_end = alignment['subject_end']
    subject_length = alignment['subject_length']
    
    # Check existing features on same chromosome and strand
    if subject in features and strand in features[subject]:
        # Sort existing features by start position
        existing_starts = sorted(features[subject][strand].keys())
        
        for prev_beg in existing_starts:
            # Skip if too far away
            distance = subject_beg - prev_beg
            if distance > 1000:
                continue
            if distance < -1000:
                break
            
            # Check each query at this start position
            prev_features = features[subject][strand][prev_beg]
            for prev_query in list(prev_features.keys()):
                prev_end = prev_features[prev_query]['subject_end']
                
                # Check for overlap with extended region
                flank_beg = max(1, subject_beg - 30)
                flank_end = min(subject_length, subject_end + 30)
                
                if overlapping(flank_beg, flank_end, prev_beg, prev_end):
                    # Merge overlapping regions
                    new_beg = min(subject_beg, prev_beg)
                    new_end = max(subject_end, prev_end)
                    
                    # Update current alignment
                    subject_beg = new_beg
                    subject_end = new_end
                    
                    # Remove merged feature
                    del prev_features[prev_query]
                    if not prev_features:
                        del features[subject][strand][prev_beg]
    
    # Store the (possibly merged) alignment
    if subject not in features:
        features[subject] = {}
    if strand not in features[subject]:
        features[subject][strand] = {}
    if subject_beg not in features[subject][strand]:
        features[subject][strand][subject_beg] = {}
    
    features[subject][strand][subject_beg][alignment['query']] = {
        'subject_end': subject_end,
        'subject_length': subject_length,
        'query_beg': alignment['query_beg'],
        'query_end': alignment['query_end'],
        'query_length': alignment['query_length']
    }

def overlapping(beg1: int, end1: int, beg2: int, end2: int) -> bool:
    """
    Check if two regions overlap.
    
    Args:
        beg1, end1: Start and end of first region
        beg2, end2: Start and end of second region
        
    Returns:
        True if regions overlap or are adjacent
    """
    return (beg1 <= beg2 <= end1 + 1) or (beg1 <= end2 + 1 <= end1) or \
           (beg2 <= beg1 <= end2 + 1) or (beg2 <= end1 + 1 <= end2)

def contained(beg1: int, end1: int, beg2: int, end2: int) -> bool:
    """
    Check if first region is contained within second.
    
    Args:
        beg1, end1: Start and end of first region
        beg2, end2: Start and end of second region
        
    Returns:
        True if first region is contained in second
    """
    return beg2 <= beg1 and end1 <= end2

def excise_candidates(features: Dict, genome_seqs: Dict, precursor_length: int, 
                      output_file: str) -> int:
    """
    Excise candidate precursor sequences and write to output.
    
    Args:
        features: Dictionary of merged features
        genome_seqs: Dictionary of genome sequences
        precursor_length: Maximum precursor length
        output_file: Path to output FASTA file
        
    Returns:
        Number of candidate precursors excised
    """
    candidate_count = 0
    short_length = precursor_length - 23
    
    try:
        with open(output_file, 'w') as out_f:
            # Sort subjects for consistent output
            sorted_subjects = sorted(features.keys())
            
            for subject in sorted_subjects:
                if subject not in genome_seqs:
                    sys.stderr.write(f"Warning: Subject {subject} not found in genome sequences\n")
                    continue
                
                genome_seq = genome_seqs[subject]
                seq_length = len(genome_seq)
                
                # Process each strand
                if '+' in features[subject]:
                    for subject_beg in sorted(features[subject]['+'].keys()):
                        for query in features[subject]['+'][subject_beg].keys():
                            feat = features[subject]['+'][subject_beg][query]
                            subject_end = feat['subject_end']
                            
                            candidate_count = excise_single_candidate(
                                subject, '+', genome_seq, seq_length,
                                subject_beg, subject_end, precursor_length,
                                short_length, candidate_count, out_f
                            )
                
                if '-' in features[subject]:
                    # For minus strand, we need to convert coordinates
                    # The input coordinates are in positive strand coordinates
                    # But for extraction on negative strand, we need to:
                    # 1. Take reverse complement of the sequence
                    # 2. Extract using a different strategy (from end backward)
                    
                    # Get the reverse complement of the sequence
                    rev_seq = revcom(genome_seq)
                    rev_length = len(rev_seq)
                    
                    for subject_beg in sorted(features[subject]['-'].keys()):
                        for query in features[subject]['-'][subject_beg].keys():
                            feat = features[subject]['-'][subject_beg][query]
                            subject_end = feat['subject_end']
                            
                            # Convert positive strand coordinates to negative strand coordinates
                            # Negative strand coordinates = L - positive_end + 1 .. L - positive_beg + 1
                            neg_beg = seq_length - subject_end + 1
                            neg_end = seq_length - subject_beg + 1
                            
                            # Now extract using negative strand coordinates on reverse complemented sequence
                            # The extraction logic should be the same as for positive strand
                            # because we're now working in the negative strand coordinate system
                            candidate_count = excise_single_candidate(
                                subject, '-', rev_seq, rev_length,
                                neg_beg, neg_end, precursor_length,
                                short_length, candidate_count, out_f,
                                # Additional parameter to indicate we need to convert coordinates back for output
                                is_negative_strand=True,
                                original_seq_length=seq_length
                            )
        
        print(f"Excised {candidate_count} candidate precursors")
        return candidate_count
        
    except Exception as e:
        sys.stderr.write(f"Error: Failed to write output file {output_file}: {e}\n")
        sys.exit(1)

def excise_single_candidate(subject: str, strand: str, seq: str, seq_length: int,
                            subject_beg: int, subject_end: int, precursor_length: int,
                            short_length: int, candidate_count: int, out_f,
                            is_negative_strand: bool = False, 
                            original_seq_length: int = None) -> int:
    """
    Excise a single candidate precursor.
    
    Args:
        subject: Subject/chromosome ID
        strand: Strand (+ or -)
        seq: Sequence to excise from
        seq_length: Length of sequence
        subject_beg: Start position of feature
        subject_end: End position of feature
        precursor_length: Maximum precursor length
        short_length: Length for short excisions
        candidate_count: Current candidate count
        out_f: Output file handle
        is_negative_strand: Whether this is a negative strand feature
        original_seq_length: Original sequence length (for negative strand coordinate conversion)
        
    Returns:
        Updated candidate count
    """
    # For longer alignments (>30 bp), excise as single precursor
    if (subject_end - subject_beg) > 30:
        excise_beg = max(1, subject_beg - 22)
        excise_end = min(seq_length, subject_end + 22)
        
        # Extend to approximately precursor_length if needed
        # This matches Perl script behavior
        current_length = excise_end - excise_beg + 1
        if current_length < precursor_length:
            # Calculate how much we need to extend
            extension_needed = precursor_length - current_length
            
            # Try to extend both sides equally
            left_ext = extension_needed // 2
            right_ext = extension_needed - left_ext
            
            excise_beg = max(1, excise_beg - left_ext)
            excise_end = min(seq_length, excise_end + right_ext)
        
        excise_lng = excise_end - excise_beg + 1
        
        # Allow some extra length beyond precursor_length
        if excise_lng <= precursor_length + 30:
            seq_sub = seq[excise_beg-1:excise_end]
            
            # For negative strand, we need to convert coordinates back to positive strand
            if is_negative_strand and original_seq_length is not None:
                # Convert negative strand coordinates back to positive strand coordinates
                pos_excise_beg = original_seq_length - excise_end + 1
                pos_excise_end = original_seq_length - excise_beg + 1
                out_f.write(f">{subject}_{candidate_count} strand:{strand} excise_beg:{pos_excise_beg} excise_end:{pos_excise_end}\n")
            else:
                out_f.write(f">{subject}_{candidate_count} strand:{strand} excise_beg:{excise_beg} excise_end:{excise_end}\n")
            
            out_f.write(f"{seq_sub}\n")
            candidate_count += 1
    else:
        # For shorter alignments, excise two precursors
        # First: centered on beginning
        excise_beg1 = max(1, subject_beg - 22)
        excise_end1 = min(seq_length, subject_beg + short_length)
        excise_lng1 = excise_end1 - excise_beg1 + 1
        
        if excise_lng1 <= precursor_length + 30:
            seq_sub1 = seq[excise_beg1-1:excise_end1]
            
            # For negative strand, convert coordinates back
            if is_negative_strand and original_seq_length is not None:
                pos_excise_beg1 = original_seq_length - excise_end1 + 1
                pos_excise_end1 = original_seq_length - excise_beg1 + 1
                out_f.write(f">{subject}_{candidate_count} strand:{strand} excise_beg:{pos_excise_beg1} excise_end:{pos_excise_end1}\n")
            else:
                out_f.write(f">{subject}_{candidate_count} strand:{strand} excise_beg:{excise_beg1} excise_end:{excise_end1}\n")
            
            out_f.write(f"{seq_sub1}\n")
            candidate_count += 1
        
        # Second: centered on end
        excise_beg2 = max(1, subject_end - short_length)
        excise_end2 = min(seq_length, subject_end + 22)
        excise_lng2 = excise_end2 - excise_beg2 + 1
        
        if excise_lng2 <= precursor_length + 30:
            seq_sub2 = seq[excise_beg2-1:excise_end2]
            
            # For negative strand, convert coordinates back
            if is_negative_strand and original_seq_length is not None:
                pos_excise_beg2 = original_seq_length - excise_end2 + 1
                pos_excise_end2 = original_seq_length - excise_beg2 + 1
                out_f.write(f">{subject}_{candidate_count} strand:{strand} excise_beg:{pos_excise_beg2} excise_end:{pos_excise_end2}\n")
            else:
                out_f.write(f">{subject}_{candidate_count} strand:{strand} excise_beg:{excise_beg2} excise_end:{excise_end2}\n")
            
            out_f.write(f"{seq_sub2}\n")
            candidate_count += 1
    
    return candidate_count

def main():
    parser = argparse.ArgumentParser(
        description='Excise potential miRNA precursor sequences from genome',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script excises potential miRNA precursor sequences from a genome
using the positions of aligned reads as guidelines. It merges nearby
alignments and extracts sequences of specified maximum length.

Examples:
  # Basic usage with default precursor length (300)
  python excise_candidate.py genome.fa alignments.bst precursors.fa
  
  # Specify precursor length
  python excise_candidate.py genome.fa alignments.bst -l 350 -o precursors.fa
        """
    )
    
    parser.add_argument('genome_fasta',
                       help='Genome FASTA file')
    
    parser.add_argument('blast_file',
                       help='BLAST-parsed alignment file')
    
    parser.add_argument('-l', '--precursor-length', type=int, default=300,
                       help='Maximum length of excised precursors (default: 300)')
    
    parser.add_argument('-o', '--output', required=True,
                       help='Output FASTA file for candidate precursors')
    
    parser.add_argument('--stats', action='store_true',
                       help='Show detailed statistics')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("EXCISING CANDIDATE miRNA PRECURSORS")
    print("=" * 60)
    print(f"Genome file: {args.genome_fasta}")
    print(f"Alignment file: {args.blast_file}")
    print(f"Maximum precursor length: {args.precursor_length}")
    print(f"Output file: {args.output}")
    print("-" * 60)
    
    # Step 1: Parse genome FASTA
    print("1. Parsing genome FASTA file...")
    genome_seqs = parse_fasta(args.genome_fasta)
    
    # Step 2: Parse BLAST alignments
    print("2. Parsing BLAST alignment file...")
    alignments = parse_blast_file(args.blast_file)
    
    # Step 3: Merge alignments into features
    print("3. Merging alignments into features...")
    features = {}
    for alignment in alignments:
        merge_features(features, alignment)
    
    # Count features
    feature_count = 0
    for subject in features.values():
        for strand_features in subject.values():
            for start_features in strand_features.values():
                feature_count += len(start_features)
    
    print(f"   Created {feature_count} features from {len(alignments)} alignments")
    
    # Step 4: Excise candidate precursors
    print("4. Excising candidate precursors...")
    candidate_count = excise_candidates(
        features, genome_seqs, args.precursor_length, args.output
    )
    
    # Statistics
    if args.stats:
        print("\n" + "=" * 60)
        print("STATISTICS")
        print("=" * 60)
        print(f"Genome sequences: {len(genome_seqs)}")
        print(f"Alignments parsed: {len(alignments)}")
        print(f"Features created: {feature_count}")
        print(f"Precursor candidates excised: {candidate_count}")
        
        # Distribution of features per chromosome
        print(f"\nFeatures per chromosome (top 10):")
        chrom_counts = {}
        for subject, strands in features.items():
            total = 0
            for strand_features in strands.values():
                for start_features in strand_features.values():
                    total += len(start_features)
            chrom_counts[subject] = total
        
        sorted_chroms = sorted(chrom_counts.items(), key=lambda x: x[1], reverse=True)[:10]
        for chrom, count in sorted_chroms:
            print(f"  {chrom}: {count} features")
    
    print("\n" + "=" * 60)
    print("EXTRACTION COMPLETE")
    print(f"Candidate precursors written to: {args.output}")
    print("=" * 60)

if __name__ == "__main__":
    main()
