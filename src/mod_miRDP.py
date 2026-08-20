#!/usr/bin/env python3
"""
mod_miRDP.py - Core miRNA prediction algorithm
"""

import sys
import os
import re
import math
from collections import defaultdict
from typing import Dict, List, Tuple, Set, Optional
import argparse
import subprocess
import tempfile

class RNAfoldResult:
    """Represents RNAfold output for a single precursor"""
    
    def __init__(self, seq_id: str, desc: str, sequence: str, structure: str, mfe: float):
        self.seq_id = seq_id
        self.desc = desc
        self.sequence = sequence
        self.structure = structure
        self.mfe = mfe
        
    def __repr__(self):
        return f"RNAfoldResult({self.seq_id}, MFE={self.mfe:.2f})"

class Alignment:
    """Represents a BLAST-parsed alignment (Perl-compatible parsing)"""
    
    def __init__(self, line: str):
        self.raw_line = line.rstrip('\n')
        parts = line.strip().split('\t')
        
        # Perl expects at least 10 fields, but we'll be flexible
        if len(parts) < 10:
            raise ValueError(f"Invalid BLAST line: {line}")
        
        # Parse based on your input format
        self.query = parts[0]
        self.query_length = int(parts[1])
        self.subject = parts[3]
        self.subject_length = int(parts[4])
        self.e_value = parts[6]
        self.identity = float(parts[7])
        self.bitscore = float(parts[8])
        self.other_info = parts[9]
        
        # Parse ranges (format: "1..24")
        self.query_beg, self.query_end = self._parse_range(parts[2])
        self.subject_beg, self.subject_end = self._parse_range(parts[5])
        
        # Determine strand
        self.strand = self._determine_strand()
        
        # Extract frequency from query ID
        self.frequency = self._extract_frequency()
    
    def _parse_range(self, range_str: str) -> Tuple[int, int]:
        """Parse range like '1..24' into (start, end)"""
        if '..' in range_str:
            beg_str, end_str = range_str.split('..')
            return int(beg_str), int(end_str)
        elif '.' in range_str:
            beg_str, end_str = range_str.split('.')
            return int(beg_str), int(end_str)
        else:
            raise ValueError(f"Invalid range format: {range_str}")
    
    def _determine_strand(self) -> str:
        """Determine strand from other_info"""
        if 'Minus' in self.other_info or '/ Minus' in self.other_info:
            return '-'
        return '+'
    
    def _extract_frequency(self) -> int:
        """Extract read frequency from query ID"""
        match = re.search(r'_x(\d+)$', self.query)
        if match:
            return int(match.group(1))
        return 1
    
    def __repr__(self):
        return f"Alignment({self.query} -> {self.subject}:{self.subject_beg}-{self.subject_end})"

class HairpinComponent:
    """Represents a component of a hairpin (mature, star, loop, flanks)"""
    
    def __init__(self):
        self.beg = 0
        self.end = 0
        self.sequence = ""
        self.structure = ""
        self.arm = ""  # "first" or "second"

class HairpinAnalysis:
    """Complete analysis of a potential miRNA hairpin"""
    
    def __init__(self, subject: str):
        self.subject = subject
        
        # Basic precursor info
        self.pri_id = ""
        self.pri_seq = ""
        self.pri_struct = ""
        self.pri_mfe = 0.0
        self.pri_beg = 1
        self.pri_end = 0
        
        # Components
        self.mature = HairpinComponent()
        self.star = HairpinComponent()
        self.loop = HairpinComponent()
        self.flank_first = HairpinComponent()
        self.flank_second = HairpinComponent()
        
        # Stem information for Drosha scoring
        self.stem_first = ""
        self.stem_second = ""
        self.stem_bp_first = 0
        self.stem_bp_second = 0
        self.stem_bp = 0
        
        # Additional info
        self.star_read_present = False
        self.total_frequency = 0
        self.score = 0.0
        self.filter_message = ""
        self.score_message = ""
        
        # Base pairing information
        self.base_pairs = {}
        
        # Read alignments for this precursor
        self.alignments = []
        self.alignment_lines = []
        self.mature_query = ""
        self.mature_strand = ""
        
        # For assembly test
        self.pre_struct = ""
        self.pre_seq = ""
    
    def calculate_base_pairs(self):
        """Calculate base pairing from structure"""
        self.base_pairs = {}
        stack = []
        
        for i, char in enumerate(self.pri_struct, 1):
            if char == '(':
                stack.append(i)
            elif char == ')':
                if stack:
                    j = stack.pop()
                    self.base_pairs[i] = j
                    self.base_pairs[j] = i
    
    def get_sequence_segment(self, beg: int, end: int, strand: str = '+') -> str:
        """Get a segment of the precursor sequence"""
        if beg < 1 or end > len(self.pri_seq) or beg > end:
            return ""
        
        seq = self.pri_seq[beg-1:end]
        if strand == '-':
            seq = self.reverse_complement(seq)
        return seq
    
    def get_structure_segment(self, beg: int, end: int, strand: str = '+') -> str:
        """Get a segment of the precursor structure"""
        if beg < 1 or end > len(self.pri_struct) or beg > end:
            return ""
        
        struct = self.pri_struct[beg-1:end]
        return struct
    
    @staticmethod
    def reverse_complement(seq: str) -> str:
        """Return reverse complement of DNA sequence"""
        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 
                'N': 'N', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n',
                'U': 'A', 'u': 'a'}
        return ''.join(comp.get(base, 'N') for base in reversed(seq))
    
    def __repr__(self):
        return f"HairpinAnalysis({self.subject}, score={self.score:.2f})"

class miRDPPredictor:
    """Main miRNA prediction class"""
    
    def __init__(self, score_threshold: float = 1.0, 
                 sensitive_mode: bool = False,
                 use_randfold: bool = False,
                 consider_drosha: bool = False,
                 debug_mode: bool = False):
        
        # Parameters
        self.score_threshold = score_threshold
        self.sensitive_mode = sensitive_mode
        self.use_randfold = use_randfold
        self.consider_drosha = consider_drosha
        self.debug_mode = debug_mode
        
        # Constants from original script
        self.nucleus_length = 7
        
        # Scoring parameters
        self.score_star = 3.9
        self.score_star_not = -1.3
        self.score_nucleus = 3.0
        self.score_nucleus_not = -0.6
        self.score_randfold = 1.6
        self.score_randfold_not = -2.2
        self.score_intercept = 0.3
        self.scores_stem = [-3.1, -2.3, -2.2, -1.6, -1.5, 0.1, 0.6, 0.8, 0.9, 0.9, 0]
        
        # Known miRNAs for conservation scoring
        self.known_nuclei = set()
        self.known_mir_info = defaultdict(list)
        
        # Data storage
        self.rnafold_results = {}
        self.alignments_by_subject = defaultdict(list)
        
        # Mathematical constant
        self.e = 2.718281828459045
        
        # For output
        self.filtered_output = False
        self.limited_output = False
        
        # Debugging
        self.skipped_subjects = []
        self.processed_count = 0
        self.passed_count = 0
    
    def parse_rnafold_output(self, rnafold_file: str):
        """Parse RNAfold output file - FIXED for your format"""
        print(f"Parsing RNAfold output: {rnafold_file}")
        
        current_id = None
        current_desc = ""
        current_seq = ""
        current_struct = ""
        current_mfe = 0.0
        
        try:
            with open(rnafold_file, 'r') as f:
                lines = f.readlines()
            
            i = 0
            while i < len(lines):
                line = lines[i].rstrip('\n')
                
                if line.startswith('>'):
                    # Save previous record
                    if current_id is not None:
                        # Apply Perl's transformations
                        current_seq = current_seq.upper().replace('U', 'T')
                        self.rnafold_results[current_id] = RNAfoldResult(
                            current_id, current_desc, current_seq, 
                            current_struct, current_mfe
                        )
                    
                    # Start new record
                    header = line[1:].strip()
                    if ' ' in header:
                        current_id, current_desc = header.split(' ', 1)
                    else:
                        current_id = header
                        current_desc = ""
                    
                    current_seq = ""
                    current_struct = ""
                    current_mfe = 0.0
                    i += 1
                    continue
                
                # Skip empty lines
                if not line.strip():
                    i += 1
                    continue
                
                # Check if line contains only sequence characters (ACGTU or N)
                if re.match(r'^[ACGTUNacgtun]+$', line):
                    # Sequence line - convert U to T as Perl does
                    current_seq += line.upper().replace('U', 'T')
                    i += 1
                else:
                    # This is the structure line which may contain MFE at the end
                    # Extract structure (dots and brackets only)
                    struct_part = re.sub(r'[^\.\(\)]', '', line)
                    current_struct = struct_part
                    
                    # Try to extract MFE from the end of the line
                    # Look for pattern like "(-62.40)"
                    mfe_match = re.search(r'\(([-\d\.]+)\)\s*$', line)
                    if mfe_match:
                        try:
                            current_mfe = float(mfe_match.group(1))
                        except ValueError:
                            current_mfe = 0.0
                    else:
                        current_mfe = 0.0
                    
                    i += 1
            
            # Save last record
            if current_id is not None:
                current_seq = current_seq.upper().replace('U', 'T')
                self.rnafold_results[current_id] = RNAfoldResult(
                    current_id, current_desc, current_seq, 
                    current_struct, current_mfe
                )
            
            print(f"  Parsed {len(self.rnafold_results)} RNAfold results")
            
        except Exception as e:
            sys.stderr.write(f"Error parsing RNAfold file: {e}\n")
            import traceback
            traceback.print_exc()
            sys.exit(1)
    
    def parse_signature_file(self, signature_file: str):
        """Parse signature file (BLAST-parsed format)"""
        print(f"Parsing signature file: {signature_file}")
        
        try:
            alignment_count = 0
            with open(signature_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:
                        continue
                    
                    try:
                        alignment = Alignment(line)
                        self.alignments_by_subject[alignment.subject].append(alignment)
                        alignment_count += 1
                    except ValueError as e:
                        if self.debug_mode:
                            print(f"  Skipping line {line_num}: {e}")
                        continue
            
            print(f"  Parsed {alignment_count} alignments for {len(self.alignments_by_subject)} subjects")
            
        except Exception as e:
            sys.stderr.write(f"Error parsing signature file: {e}\n")
            sys.exit(1)
    
    def parse_known_mirnas(self, mirna_file: str):
        """Parse known miRNA sequences for conservation scoring"""
        print(f"Parsing known miRNAs: {mirna_file}")
        
        try:
            with open(mirna_file, 'r') as f:
                current_id = ""
                current_seq = ""
                
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    if line.startswith('>'):
                        # Process previous sequence
                        if current_id and current_seq:
                            if len(current_seq) >= self.nucleus_length + 1:
                                nucleus = current_seq[1:self.nucleus_length + 1]
                                nucleus = nucleus.upper().replace('T', 'U')
                                self.known_nuclei.add(nucleus)
                                self.known_mir_info[nucleus].append(current_id)
                        
                        # Start new sequence
                        current_id = line[1:].split()[0]
                        current_seq = ""
                    else:
                        current_seq += line.upper()
                
                # Process last sequence
                if current_id and current_seq:
                    if len(current_seq) >= self.nucleus_length + 1:
                        nucleus = current_seq[1:self.nucleus_length + 1]
                        nucleus = nucleus.upper().replace('T', 'U')
                        self.known_nuclei.add(nucleus)
                        self.known_mir_info[nucleus].append(current_id)
            
            print(f"  Loaded {len(self.known_nuclei)} unique seed sequences")
            
        except Exception as e:
            sys.stderr.write(f"Error parsing known miRNA file: {e}\n")
            sys.exit(1)
    
    def run_randfold(self, sequence: str) -> float:
        """Run Randfold to get p-value for sequence"""
        try:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp:
                tmp.write(f">pri_seq\n{sequence}\n")
                tmp_name = tmp.name
            
            cmd = ['randfold', '-s', tmp_name, '999']
            result = subprocess.run(cmd, capture_output=True, text=True, shell=False)
            
            os.unlink(tmp_name)
            
            if result.returncode == 0:
                lines = result.stdout.strip().split('\n')
                if len(lines) >= 2:
                    fields = lines[1].split()
                    if len(fields) >= 3:
                        return float(fields[2])
            
            return 1.0
            
        except Exception as e:
            if self.debug_mode:
                sys.stderr.write(f"Warning: Randfold failed: {e}\n")
            return 1.0
    
    def analyze_hairpin(self, subject: str) -> Optional[HairpinAnalysis]:
        """Analyze a single potential hairpin"""
        if self.debug_mode:
            print(f"DEBUG: Analyzing hairpin: {subject}")
        
        # Try exact match first
        if subject in self.rnafold_results:
            rnafold = self.rnafold_results[subject]
            actual_subject = subject
        else:
            # Try to find subject without suffix
            found = False
            actual_subject = None
            
            # Try removing numeric suffix
            base_id = re.sub(r'_\d+$', '', subject)
            if base_id in self.rnafold_results:
                rnafold = self.rnafold_results[base_id]
                actual_subject = base_id
                found = True
            
            # Try common suffixes
            if not found:
                for suffix in ['_1', '_2', '_3', '_4', '_5']:
                    test_id = subject + suffix
                    if test_id in self.rnafold_results:
                        rnafold = self.rnafold_results[test_id]
                        actual_subject = test_id
                        found = True
                        break
            
            if not found:
                if self.debug_mode:
                    print(f"  No RNAfold result found for {subject}")
                    self.skipped_subjects.append(subject)
                return None
        
        # Check if we have structure
        if not rnafold.structure:
            if self.debug_mode:
                print(f"  No structure for {actual_subject}")
            return None
        
        # Create HairpinAnalysis object
        hairpin = HairpinAnalysis(actual_subject)
        hairpin.pri_id = actual_subject
        hairpin.pri_seq = rnafold.sequence
        hairpin.pri_struct = rnafold.structure
        hairpin.pri_mfe = rnafold.mfe
        hairpin.pri_end = len(rnafold.sequence)
        
        # Get alignments for this subject
        alignments = []
        if subject in self.alignments_by_subject:
            alignments = self.alignments_by_subject[subject]
        elif actual_subject in self.alignments_by_subject:
            alignments = self.alignments_by_subject[actual_subject]
        else:
            # Try other variations
            for key, aligns in self.alignments_by_subject.items():
                if key == subject or key == actual_subject or key.startswith(subject + '_'):
                    alignments = aligns
                    break
        
        hairpin.alignments = alignments
        hairpin.alignment_lines = [align.raw_line for align in alignments]
        
        if not hairpin.alignments:
            if self.debug_mode:
                print(f"  No alignments for {actual_subject}")
            return None
        
        # Calculate base pairing
        hairpin.calculate_base_pairs()
        
        # Find mature sequence (read with highest frequency)
        hairpin.alignments.sort(key=lambda x: (x.frequency, x.query), reverse=True)
        mature_alignment = hairpin.alignments[0]
        
        hairpin.mature.beg = mature_alignment.subject_beg
        hairpin.mature.end = mature_alignment.subject_end
        hairpin.mature.sequence = hairpin.get_sequence_segment(
            mature_alignment.subject_beg, mature_alignment.subject_end, mature_alignment.strand
        )
        hairpin.mature.structure = hairpin.get_structure_segment(
            mature_alignment.subject_beg, mature_alignment.subject_end, mature_alignment.strand
        )
        
        hairpin.mature_query = mature_alignment.query
        hairpin.mature_strand = mature_alignment.strand
        
        # Determine mature arm
        hairpin.mature.arm = self._determine_arm(hairpin, 
            mature_alignment.subject_beg, mature_alignment.subject_end, mature_alignment.strand)
        
        if not hairpin.mature.arm:
            if self.debug_mode:
                print(f"  Could not determine mature arm for {actual_subject}")
            return None
        
        # Find star sequence
        self._find_star_sequence(hairpin)
        
        # Find loop sequence
        self._find_loop_sequence(hairpin)
        
        # Find flanking sequences
        self._find_flanking_sequences(hairpin)
        
        # Assemble precursor structure and sequence
        self._assemble_precursor(hairpin)
        
        if self.debug_mode:
            print(f"  Successfully analyzed {actual_subject}")
        
        return hairpin
    
    def _determine_arm(self, hairpin: HairpinAnalysis, beg: int, end: int, strand: str) -> str:
        """Determine if sequence is in first (5') or second (3') arm"""
        if strand == '-':
            return ""
        
        struct = hairpin.get_structure_segment(beg, end, strand)
        
        if not struct:
            return ""
        
        if re.match(r'^[\.\(]+$', struct) and '(' in struct:
            return "first"
        elif re.match(r'^[\.\)]+$', struct) and ')' in struct:
            return "second"
        
        return ""
    
    def _find_star_sequence(self, hairpin: HairpinAnalysis):
        """Find star sequence based on base pairing with mature"""
        if not hairpin.mature.arm:
            return
        
        mature_beg = hairpin.mature.beg
        mature_end = hairpin.mature.end
        
        # Consider 2nt 3' overhang
        mature_end_for_star = max(mature_beg, mature_end - 2)
        mature_lng = mature_end_for_star - mature_beg + 1
        
        # Find star beginning
        star_beg = None
        offset_beg = 0
        while not star_beg and offset_beg < mature_lng:
            pos = mature_end_for_star - offset_beg
            if pos in hairpin.base_pairs:
                star_beg_pos = hairpin.base_pairs[pos]
                star_beg = star_beg_pos - offset_beg
                break
            offset_beg += 1
        
        # Find star end
        star_end = None
        offset_end = 0
        while not star_end and offset_end < mature_lng:
            pos = mature_beg + offset_end
            if pos in hairpin.base_pairs:
                star_end_pos = hairpin.base_pairs[pos]
                star_end = star_end_pos + offset_end + 2
                break
            offset_end += 1
        
        if star_beg and star_end and star_beg <= star_end:
            star_end = min(star_end, hairpin.pri_end)
            
            hairpin.star.beg = star_beg
            hairpin.star.end = star_end
            hairpin.star.sequence = hairpin.get_sequence_segment(star_beg, star_end, '+')
            hairpin.star.structure = hairpin.get_structure_segment(star_beg, star_end, '+')
            
            struct = hairpin.star.structure
            if struct and re.match(r'^[\.\(]+$', struct) and '(' in struct:
                hairpin.star.arm = "first"
            elif struct and re.match(r'^[\.\)]+$', struct) and ')' in struct:
                hairpin.star.arm = "second"
    
    def _find_loop_sequence(self, hairpin: HairpinAnalysis):
        """Find loop sequence between mature and star"""
        if not hairpin.mature.arm or not hairpin.star.arm:
            return
        
        loop_beg = 0
        loop_end = 0
        
        if hairpin.mature.arm == "first":
            loop_beg = hairpin.mature.end + 1
        else:
            loop_end = hairpin.mature.beg - 1
        
        if hairpin.star.arm == "first":
            loop_beg = hairpin.star.end + 1
        else:
            loop_end = hairpin.star.beg - 1
        
        if (loop_beg > 0 and loop_beg <= hairpin.pri_end and 
            loop_end > 0 and loop_end <= hairpin.pri_end and 
            loop_beg <= loop_end):
            
            hairpin.loop.beg = loop_beg
            hairpin.loop.end = loop_end
            hairpin.loop.sequence = hairpin.get_sequence_segment(loop_beg, loop_end, '+')
            hairpin.loop.structure = hairpin.get_structure_segment(loop_beg, loop_end, '+')
    
    def _find_flanking_sequences(self, hairpin: HairpinAnalysis):
        """Find flanking sequences (lower stems)"""
        if not hairpin.mature.arm or not hairpin.star.arm:
            return
        
        flank_first_end = 0
        flank_second_beg = 0
        
        if hairpin.mature.arm == "first":
            flank_first_end = hairpin.mature.beg - 1
        else:
            flank_second_beg = hairpin.mature.end + 1
        
        if hairpin.star.arm == "first":
            flank_first_end = hairpin.star.beg - 1
        else:
            flank_second_beg = hairpin.star.end + 1
        
        if flank_first_end > 0:
            hairpin.flank_first.end = flank_first_end
            hairpin.flank_first.sequence = hairpin.get_sequence_segment(1, flank_first_end, '+')
            hairpin.flank_first.structure = hairpin.get_structure_segment(1, flank_first_end, '+')
        
        if flank_second_beg > 0 and flank_second_beg <= hairpin.pri_end:
            hairpin.flank_second.beg = flank_second_beg
            hairpin.flank_second.sequence = hairpin.get_sequence_segment(
                flank_second_beg, hairpin.pri_end, '+')
            hairpin.flank_second.structure = hairpin.get_structure_segment(
                flank_second_beg, hairpin.pri_end, '+')
        
        # Score stems for Drosha recognition
        if self.consider_drosha:
            self._score_stems_drosha(hairpin)
    
    def _score_stems_drosha(self, hairpin: HairpinAnalysis):
        """Score lower stems for Drosha recognition"""
        if not hairpin.flank_first.structure or not hairpin.flank_second.structure:
            return
        
        stem_first = hairpin.flank_first.structure[-10:] if len(hairpin.flank_first.structure) >= 10 else hairpin.flank_first.structure
        stem_second = hairpin.flank_second.structure[:10] if len(hairpin.flank_second.structure) >= 10 else hairpin.flank_second.structure
        
        hairpin.stem_first = stem_first
        hairpin.stem_second = stem_second
        
        hairpin.stem_bp_first = stem_first.count('(')
        hairpin.stem_bp_second = stem_second.count(')')
        
        hairpin.stem_bp = min(hairpin.stem_bp_first, hairpin.stem_bp_second)
    
    def _assemble_precursor(self, hairpin: HairpinAnalysis):
        """Assemble precursor from components"""
        if not hairpin.mature.arm or not hairpin.star.arm:
            return
        
        if hairpin.mature.arm == "first":
            hairpin.pre_seq = hairpin.mature.sequence + hairpin.loop.sequence + hairpin.star.sequence
            hairpin.pre_struct = hairpin.mature.structure + hairpin.loop.structure + hairpin.star.structure
        else:
            hairpin.pre_seq = hairpin.star.sequence + hairpin.loop.sequence + hairpin.mature.sequence
            hairpin.pre_struct = hairpin.star.structure + hairpin.loop.structure + hairpin.mature.structure
    
    def filter_initial(self, hairpin: HairpinAnalysis) -> bool:
        """Apply initial filtering criteria"""
        hairpin.filter_message = ""
        
        if not self._pass_filtering_structure(hairpin):
            if hairpin.filter_message:
                hairpin.filter_message = "structure problem\n" + hairpin.filter_message
            else:
                hairpin.filter_message = "structure problem\n"
            return False
        
        if not self._pass_filtering_signature(hairpin):
            if hairpin.filter_message:
                hairpin.filter_message = "signature problem\n" + hairpin.filter_message
            else:
                hairpin.filter_message = "signature problem\n"
            return False
        
        return True
    
    def _pass_filtering_structure(self, hairpin: HairpinAnalysis) -> bool:
        """Test structure"""
        if not self._test_components(hairpin):
            return False
        
        if not self._check_no_bifurcations(hairpin):
            return False
        
        if self._count_base_pairs_duplex(hairpin) < 14:
            hairpin.filter_message += "too few pairings in duplex\n"
            return False
        
        mature_len = len(hairpin.mature.sequence) if hairpin.mature.sequence else 0
        star_len = len(hairpin.star.sequence) if hairpin.star.sequence else 0
        if abs(mature_len - star_len) >= 6:
            hairpin.filter_message += "too big difference between mature and star length\n"
            return False
        
        return True
    
    def _test_components(self, hairpin: HairpinAnalysis) -> bool:
        """Test if all components are present"""
        if not hairpin.mature.structure:
            hairpin.filter_message += "no mature\n"
            return False
        
        if not hairpin.star.structure:
            hairpin.filter_message += "no star\n"
            return False
        
        if not hairpin.loop.structure:
            hairpin.filter_message += "no loop\n"
            return False
        
        if not hairpin.flank_first.structure:
            hairpin.filter_message += "no flanks\n"
            return False
        
        if not hairpin.flank_second.structure:
            hairpin.filter_message += "no flanks\n"
            return False
        
        return True
    
    def _check_no_bifurcations(self, hairpin: HairpinAnalysis) -> bool:
        """Check for bifurcations"""
        mature_struc1 = ""
        mature_struc2 = ""
        
        if hairpin.mature.arm == "first":
            mature_struc1 = hairpin.mature.structure + hairpin.loop.structure[:5]
            mature_struc2 = hairpin.loop.structure[-5:] + hairpin.star.structure
        else:
            mature_struc1 = hairpin.star.structure + hairpin.loop.structure[:5]
            mature_struc2 = hairpin.loop.structure[-5:] + hairpin.mature.structure
        
        if len(hairpin.loop.structure) >= 20:
            if not ((re.match(r'^[\.\(]+$', mature_struc1) or re.match(r'^[\.\)]+$', mature_struc1)) and
                    (re.match(r'^[\.\(]+$', mature_struc2) or re.match(r'^[\.\)]+$', mature_struc2))):
                hairpin.filter_message += "bifurcation in precursor\n"
                return False
        else:
            if not ((re.match(r'^[\.\(]+$', hairpin.mature.structure) or re.match(r'^[\.\)]+$', hairpin.mature.structure)) and
                    (re.match(r'^[\.\(]+$', hairpin.star.structure) or re.match(r'^[\.\)]+$', hairpin.star.structure))):
                hairpin.filter_message += "bifurcation in precursor\n"
                return False
        
        return True
    
    def _count_base_pairs_duplex(self, hairpin: HairpinAnalysis) -> int:
        """Count base pairs in duplex"""
        mature_struct = hairpin.mature.structure
        if not mature_struct:
            return 0
        return mature_struct.count('(') + mature_struct.count(')')
    
    def _pass_filtering_signature(self, hairpin: HairpinAnalysis) -> bool:
        """Check read distribution"""
        consistent = 0
        inconsistent = 0
        star_perfect = 0
        star_fuzzy = 0
        
        sorted_alignments = sorted(hairpin.alignments, key=lambda x: x.subject_beg)
        
        for alignment in sorted_alignments:
            freq = alignment.frequency
            product = self._classify_read(alignment, hairpin)
            
            if product:
                consistent += freq
                if product == "star" and self._check_star_overhang(alignment, hairpin):
                    star_perfect += freq
                elif product == "star":
                    star_fuzzy += freq
            else:
                inconsistent += freq
        
        total = consistent + inconsistent
        hairpin.total_frequency = total
        
        if total == 0:
            hairpin.filter_message += "read frequency too low\n"
            return False
        
        if star_perfect > star_fuzzy:
            hairpin.star_read_present = True
        
        if inconsistent > 0:
            inconsistent_fraction = inconsistent / total
            if inconsistent_fraction > 0.2:
                hairpin.filter_message += f"inconsistent\t{inconsistent}\nconsistent\t{consistent}\n"
                return False
        
        return True
    
    def _classify_read(self, alignment: Alignment, hairpin: HairpinAnalysis) -> Optional[str]:
        """Classify a read"""
        if alignment.strand == '-':
            return None
        
        beg = alignment.subject_beg
        end = alignment.subject_end
        
        fuzz_beg = 2
        fuzz_end = 5
        
        def contained(beg1, end1, beg2, end2):
            return beg2 <= beg1 and end1 <= end2
        
        mature_beg = hairpin.mature.beg
        mature_end = hairpin.mature.end
        if contained(beg, end, mature_beg - fuzz_beg, mature_end + fuzz_end):
            return "mature"
        
        star_beg = hairpin.star.beg
        star_end = hairpin.star.end
        if star_beg > 0 and star_end > 0:
            if contained(beg, end, star_beg - fuzz_beg, star_end + fuzz_end):
                return "star"
        
        loop_beg = hairpin.loop.beg
        loop_end = hairpin.loop.end
        if loop_beg > 0 and loop_end > 0:
            if contained(beg, end, loop_beg - fuzz_beg, loop_end + fuzz_end):
                return "loop"
        
        return None
    
    def _check_star_overhang(self, alignment: Alignment, hairpin: HairpinAnalysis) -> bool:
        """Check star overhang"""
        beg = alignment.subject_beg
        offset = beg - hairpin.star.beg
        return offset in [-1, 0, 1]
    
    def calculate_score(self, hairpin: HairpinAnalysis) -> float:
        """Calculate comprehensive score"""
        score = 0.0
        score_messages = []
        
        # MFE score
        score_mfe = self._score_mfe(hairpin.pri_mfe)
        score += score_mfe
        round_mfe = round(score_mfe * 10) / 10
        score_messages.append(f"score_mfe\t{round_mfe:.1f}")
        
        # Read frequency score
        score_freq = self._score_frequency(hairpin.total_frequency, hairpin.star_read_present)
        score += score_freq
        round_freq = round(score_freq)
        score_messages.append(f"score_freq\t{round_freq}")
        
        # Star evidence score
        if hairpin.star_read_present:
            score += self.score_star
            score_messages.append(f"score_star\t{self.score_star}")
        else:
            score += self.score_star_not
            score_messages.append(f"score_star\t{self.score_star_not}")
        
        # Conservation score
        if self.known_nuclei:
            if self._check_conservation(hairpin):
                score += self.score_nucleus
                score_messages.append(f"score_nucleus\t{self.score_nucleus}")
                if len(hairpin.mature.sequence) >= self.nucleus_length + 1:
                    nucleus = hairpin.mature.sequence[1:self.nucleus_length + 1]
                    nucleus = nucleus.upper().replace('T', 'U')
                    if nucleus in self.known_mir_info:
                        mir_ids = "\t".join(self.known_mir_info[nucleus])
                        score_messages.append(f"{mir_ids}")
            else:
                score += self.score_nucleus_not
                score_messages.append(f"score_nucleus\t{self.score_nucleus_not}")
        
        # Stem score
        if self.consider_drosha:
            stem_score = self._score_stems(hairpin)
            score += stem_score
            score_messages.append(f"score_stem\t{stem_score}")
        
        # Add intercept
        score += self.score_intercept
        
        # Randfold score
        if self.use_randfold:
            if (score + self.score_randfold >= self.score_threshold or 
                score + self.score_randfold_not <= self.score_threshold):
                
                p_value = self.run_randfold(hairpin.pri_seq)
                if p_value <= 0.05:
                    score += self.score_randfold
                    score_messages.append(f"score_randfold\t{self.score_randfold}")
                else:
                    score += self.score_randfold_not
                    score_messages.append(f"score_randfold\t{self.score_randfold_not}")
        
        # Final score
        final_score = round(score * 10) / 10
        score_messages.append(f"score\t{final_score:.1f}")
        
        hairpin.score = final_score
        hairpin.score_message = "\n".join(score_messages)
        
        return final_score
    
    def _score_mfe(self, mfe: float) -> float:
        """Score based on minimum free energy"""
        if mfe >= 0:
            return -10.0
        
        mfe_adj = max(1.0, -mfe)
        
        prob_test = self._prob_gumbel(mfe_adj, 5.5, 32)
        prob_background = self._prob_gumbel(mfe_adj, 4.8, 23)
        
        if prob_background == 0:
            return 5.4
        
        odds = prob_test / prob_background
        if odds <= 0:
            return -10.0
        
        return math.log(odds)
    
    def _prob_gumbel(self, x: float, scale: float, location: float) -> float:
        """Discretized Gumbel distribution"""
        bound_lower = x - 0.5
        bound_upper = x + 0.5
        
        cdf_lower = self._cdf_gumbel(bound_lower, scale, location)
        cdf_upper = self._cdf_gumbel(bound_upper, scale, location)
        
        return cdf_upper - cdf_lower
    
    def _cdf_gumbel(self, x: float, scale: float, location: float) -> float:
        """CDF of Gumbel distribution"""
        if scale <= 0:
            return 0.0
        return self.e ** (-(self.e ** (-(x - location) / scale)))
    
    def _score_frequency(self, frequency: int, star_read_present: bool) -> float:
        """Score based on read frequency"""
        if frequency <= 0:
            return -10.0
        
        parameter_test = 0.999
        parameter_control = 0.6
        
        intercept = math.log((1 - parameter_test) / (1 - parameter_control))
        slope = math.log(parameter_test / parameter_control)
        log_odds = slope * frequency + intercept
        
        if not self.sensitive_mode and not star_read_present:
            log_odds = min(log_odds, 0.0)
        
        return log_odds
    
    def _check_conservation(self, hairpin: HairpinAnalysis) -> bool:
        """Check conservation"""
        if len(hairpin.mature.sequence) < self.nucleus_length + 1:
            return False
        
        nucleus = hairpin.mature.sequence[1:self.nucleus_length + 1]
        nucleus = nucleus.upper().replace('T', 'U')
        return nucleus in self.known_nuclei
    
    def _score_stems(self, hairpin: HairpinAnalysis) -> float:
        """Score stems"""
        if hairpin.stem_bp < len(self.scores_stem):
            return self.scores_stem[hairpin.stem_bp]
        return 0.0
    
    def process_all(self) -> List[HairpinAnalysis]:
        """Process all potential precursors"""
        if self.debug_mode:
            print(f"\nDEBUG: Processing precursors")
        
        results = []
        processed_subjects = set()
        
        # Process subjects from signature file
        all_subjects = list(self.alignments_by_subject.keys())
        
        for i, subject in enumerate(all_subjects, 1):
            if subject in processed_subjects:
                continue
            
            if self.debug_mode and i % 100 == 0:
                print(f"  Processed {i}/{len(all_subjects)} subjects...")
            
            hairpin = self.analyze_hairpin(subject)
            self.processed_count += 1
            
            if not hairpin:
                continue
            
            if not self.filter_initial(hairpin):
                if self.filtered_output:
                    results.append(hairpin)
                continue
            
            score = self.calculate_score(hairpin)
            
            if score >= self.score_threshold:
                if not self.filtered_output:
                    self.passed_count += 1
                    results.append(hairpin)
            elif self.filtered_output:
                results.append(hairpin)
            
            processed_subjects.add(subject)
        
        if self.debug_mode:
            print(f"\nDEBUG: Processed {self.processed_count} subjects")
            print(f"DEBUG: {self.passed_count} passed threshold")
        
        return results
    
    def write_output(self, results: List[HairpinAnalysis], output_file: str):
        """Write results to output file"""
        try:
            with open(output_file, 'w') as f:
                for hairpin in results:
                    if self.limited_output:
                        f.write(f"{hairpin.pri_id}\n")
                        continue
                    
                    if hairpin.filter_message:
                        f.write(hairpin.filter_message)
                    
                    if hairpin.score_message:
                        f.write(hairpin.score_message)
                        f.write("\n")
                    
                    # Write component information
                    if hairpin.flank_first.end > 0:
                        f.write(f"flank_first_end  \t{hairpin.flank_first.end}\n")
                        f.write(f"flank_first_seq  \t{hairpin.flank_first.sequence}\n")
                        f.write(f"flank_first_struct  \t{hairpin.flank_first.structure}\n")
                    
                    if hairpin.flank_second.beg > 0:
                        f.write(f"flank_second_beg  \t{hairpin.flank_second.beg}\n")
                        f.write(f"flank_second_seq  \t{hairpin.flank_second.sequence}\n")
                        f.write(f"flank_second_struct  \t{hairpin.flank_second.structure}\n")
                    
                    f.write(f"freq  \t{hairpin.total_frequency}\n")
                    
                    if hairpin.loop.beg > 0:
                        f.write(f"loop_beg  \t{hairpin.loop.beg}\n")
                        f.write(f"loop_end  \t{hairpin.loop.end}\n")
                        f.write(f"loop_seq  \t{hairpin.loop.sequence}\n")
                        f.write(f"loop_struct  \t{hairpin.loop.structure}\n")
                    
                    if hairpin.mature.arm:
                        f.write(f"mature_arm  \t{hairpin.mature.arm}\n")
                        f.write(f"mature_beg  \t{hairpin.mature.beg}\n")
                        f.write(f"mature_end  \t{hairpin.mature.end}\n")
                        f.write(f"mature_query  \t{hairpin.mature_query}\n")
                        f.write(f"mature_seq  \t{hairpin.mature.sequence}\n")
                        f.write(f"mature_strand  \t{hairpin.mature_strand}\n")
                        f.write(f"mature_struct  \t{hairpin.mature.structure}\n")
                    
                    if hairpin.pre_seq:
                        f.write(f"pre_seq  \t{hairpin.pre_seq}\n")
                        f.write(f"pre_struct  \t{hairpin.pre_struct}\n")
                    
                    f.write(f"pri_beg  \t{hairpin.pri_beg}\n")
                    f.write(f"pri_end  \t{hairpin.pri_end}\n")
                    f.write(f"pri_id  \t{hairpin.pri_id}\n")
                    f.write(f"pri_mfe  \t{hairpin.pri_mfe:.2f}\n")
                    f.write(f"pri_seq  \t{hairpin.pri_seq}\n")
                    f.write(f"pri_struct  \t{hairpin.pri_struct}\n")
                    
                    if hairpin.star.arm:
                        f.write(f"star_arm  \t{hairpin.star.arm}\n")
                        f.write(f"star_beg  \t{hairpin.star.beg}\n")
                        f.write(f"star_end  \t{hairpin.star.end}\n")
                        f.write(f"star_seq  \t{hairpin.star.sequence}\n")
                        f.write(f"star_struct  \t{hairpin.star.structure}\n")
                    
                    if self.consider_drosha and hairpin.stem_first:
                        f.write(f"stem_first  \t{hairpin.stem_first}\n")
                        f.write(f"stem_second  \t{hairpin.stem_second}\n")
                        f.write(f"stem_bp_first  \t{hairpin.stem_bp_first}\n")
                        f.write(f"stem_bp_second  \t{hairpin.stem_bp_second}\n")
                        f.write(f"stem_bp  \t{hairpin.stem_bp}\n")
                    
                    for line in hairpin.alignment_lines:
                        f.write(f"{line}\n")
                    
                    f.write("\n")
            
            print(f"Results written to: {output_file}")
            
        except Exception as e:
            sys.stderr.write(f"Error writing output file: {e}\n")
            sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='Core miRNA prediction algorithm (Perl-compatible)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python mod_miRDP.py signatures.bst structures.txt -o predictions.txt
  
  # With debug mode
  python mod_miRDP.py signatures.bst structures.txt -o predictions.txt --debug
        """
    )
    
    parser.add_argument('signature_file',
                       help='BLAST-parsed signature file')
    
    parser.add_argument('structure_file',
                       help='RNAfold output file')
    
    parser.add_argument('-o', '--output', required=True,
                       help='Output file for predictions')
    
    parser.add_argument('-s', '--known-mirnas',
                       help='FASTA file with known miRNA sequences for conservation')
    
    parser.add_argument('-t', '--filtered-output', action='store_true',
                       help='Print precursors that do NOT pass the threshold')
    
    parser.add_argument('-u', '--limited-output', action='store_true',
                       help='Limited output (only precursor IDs)')
    
    parser.add_argument('-v', '--threshold', type=float, default=1.0,
                       help='Score threshold (default: 1.0)')
    
    parser.add_argument('-x', '--sensitive', action='store_true',
                       help='Sensitive mode for Sanger sequences')
    
    parser.add_argument('-y', '--randfold', action='store_true',
                       help='Use Randfold for additional scoring')
    
    parser.add_argument('-z', '--drosha', action='store_true',
                       help='Consider Drosha processing in scoring')
    
    parser.add_argument('--debug', action='store_true',
                       help='Enable debug output')
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("miRDP CORE ALGORITHM - PERL-COMPATIBLE VERSION")
    print("=" * 70)
    print(f"Signature file: {args.signature_file}")
    print(f"Structure file: {args.structure_file}")
    print(f"Output file:    {args.output}")
    print(f"Parameters:")
    print(f"  Score threshold: {args.threshold}")
    print(f"  Sensitive mode:  {'Yes' if args.sensitive else 'No'}")
    print(f"  Use Randfold:    {'Yes' if args.randfold else 'No'}")
    print(f"  Consider Drosha: {'Yes' if args.drosha else 'No'}")
    print(f"  Filtered output: {'Yes' if args.filtered_output else 'No'}")
    print(f"  Limited output:  {'Yes' if args.limited_output else 'No'}")
    print(f"  Debug mode:      {'Yes' if args.debug else 'No'}")
    if args.known_mirnas:
        print(f"  Known miRNAs:    {args.known_mirnas}")
    print("-" * 70)
    
    # Initialize predictor
    predictor = miRDPPredictor(
        score_threshold=args.threshold,
        sensitive_mode=args.sensitive,
        use_randfold=args.randfold,
        consider_drosha=args.drosha,
        debug_mode=args.debug
    )
    
    # Set output options
    predictor.filtered_output = args.filtered_output
    predictor.limited_output = args.limited_output
    
    # Parse known miRNAs if provided
    if args.known_mirnas:
        predictor.parse_known_mirnas(args.known_mirnas)
    
    # Parse input files
    predictor.parse_rnafold_output(args.structure_file)
    predictor.parse_signature_file(args.signature_file)
    
    # Process all precursors
    results = predictor.process_all()
    
    # Write output
    predictor.write_output(results, args.output)
    
    print("\n" + "=" * 70)
    print("PREDICTION COMPLETE")
    print("=" * 70)
    print(f"Total precursors analyzed: {predictor.processed_count}")
    if args.filtered_output:
        print(f"Precursors filtered out: {len(results)}")
    else:
        print(f"Precursors passed threshold: {len(results)}")
    print("=" * 70)

if __name__ == "__main__":
    main()
