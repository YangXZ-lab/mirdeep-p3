#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
Remove redundant miRNA predictions and apply plant-specific criteria.

This script processes miRDeep prediction results to remove identical predictions
and apply plant-specific filtering criteria.
"""

import sys
import os
import re
from collections import defaultdict
from typing import Dict, List, Tuple, Set
import argparse


def parse_chromosome_lengths(chr_length_file: str) -> Dict[str, int]:
    """Parse chromosome length file.
    
    Args:
        chr_length_file: Path to chromosome length file
        
    Returns:
        Dictionary mapping chromosome ID to length
    """
    chr_lengths = {}
    
    try:
        with open(chr_length_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) >= 2:
                    chr_lengths[parts[0]] = int(parts[1])
        
        print(f"Parsed {len(chr_lengths)} chromosomes")
        return chr_lengths
        
    except Exception as e:
        sys.stderr.write(f"Error parsing chromosome length file {chr_length_file}: {e}\n")
        sys.exit(1)


def parse_precursors(precursor_file: str) -> Dict[str, Dict]:
    """Parse precursor FASTA file and extract positional information.
    
    Args:
        precursor_file: Path to precursor FASTA file
        
    Returns:
        Dictionary with precursor information
    """
    precursors = {}
    
    try:
        with open(precursor_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('>'):
                    header = line[1:]
                    parts = header.split()
                    
                    if len(parts) < 4:
                        sys.stderr.write(f"Warning: Malformed precursor header: {line}\n")
                        continue
                    
                    precursor_id = parts[0]
                    
                    strand_part = parts[1]
                    strand = strand_part.split(':')[1] if strand_part.startswith('strand:') else '+'
                    
                    excise_beg = 0
                    excise_end = 0
                    
                    for part in parts[2:]:
                        if part.startswith('excise_beg:'):
                            try:
                                excise_beg = int(part.split(':')[1])
                            except (ValueError, IndexError):
                                excise_beg = 0
                        elif part.startswith('excise_end:'):
                            try:
                                excise_end = int(part.split(':')[1])
                            except (ValueError, IndexError):
                                excise_end = 0
                    
                    base_id = '_'.join(precursor_id.split('_')[:-1]) if '_' in precursor_id else precursor_id
                    
                    precursors[precursor_id] = {
                        'strand': strand,
                        'excise_beg': excise_beg,
                        'excise_end': excise_end,
                        'chromosome': base_id
                    }
        
        print(f"Parsed {len(precursors)} precursors")
        return precursors
        
    except Exception as e:
        sys.stderr.write(f"Error parsing precursor file {precursor_file}: {e}\n")
        sys.exit(1)


def parse_miRDeep_result(miRDeep_file: str) -> List[Dict]:
    """Parse miRDeep prediction results.
    
    Args:
        miRDeep_file: Path to miRDeep results file
        
    Returns:
        List of miRNA prediction dictionaries
    """
    predictions = []
    
    try:
        print(f"Parsing miRDeep file: {miRDeep_file}")
        
        with open(miRDeep_file, 'r') as f:
            content = f.read()
        
        blocks = []
        current_block = []
        
        lines = content.strip().split('\n')
        for line in lines:
            line = line.strip()
            
            if not line and current_block:
                blocks.append('\n'.join(current_block))
                current_block = []
            elif line:
                current_block.append(line)
        
        if current_block:
            blocks.append('\n'.join(current_block))
        
        print(f"Found {len(blocks)} blocks in miRDeep file")
        
        for block in blocks:
            if not block.strip():
                continue
            
            pred = {}
            alignment_lines = []
            
            for line in block.strip().split('\n'):
                line = line.strip()
                if not line:
                    continue
                
                if re.search(r'_x\d+\s+\d+\s+\d+\.\.\d+\s+', line):
                    alignment_lines.append(line)
                elif '\t' in line:
                    parts = line.split('\t', 1)
                    if len(parts) == 2:
                        key, value = parts
                        if not re.search(r'_x\d+$', key):
                            pred[key.strip()] = value.strip()
                elif '  ' in line:
                    parts = line.split(None, 1)
                    if len(parts) == 2:
                        key, value = parts
                        if not re.search(r'_x\d+$', key):
                            pred[key.strip()] = value
            
            if alignment_lines:
                pred['alignment_lines'] = alignment_lines
            
            if 'pri_id' in pred or 'mature_query' in pred:
                predictions.append(pred)
        
        print(f"Parsed {len(predictions)} miRNA predictions")
        return predictions
        
    except Exception as e:
        sys.stderr.write(f"Error parsing miRDeep result file {miRDeep_file}: {e}\n")
        sys.exit(1)


def get_total_reads(total_reads_file: str) -> int:
    """Read total reads count from file.
    
    Args:
        total_reads_file: Path to total reads file
        
    Returns:
        Total reads count
    """
    try:
        with open(total_reads_file, 'r') as f:
            content = f.read().strip()
            if content.endswith('%'):
                content = content[:-1]
            return int(float(content))
    except Exception as e:
        sys.stderr.write(f"Error reading total reads file {total_reads_file}: {e}\n")
        sys.exit(1)


def _find_precursor_info(pri_id: str, precursors: Dict[str, Dict]) -> Dict:
    """Helper function to find precursor information."""
    if pri_id in precursors:
        return precursors[pri_id]
    
    for key in precursors.keys():
        if key.startswith(pri_id):
            return precursors[key]
    
    for key in precursors.keys():
        base_key = '_'.join(key.split('_')[:-1]) if '_' in key else key
        if base_key == pri_id:
            return precursors[key]
    
    return None


def _validate_position(pos: int, pos_name: str, pri_id: str) -> bool:
    """Validate position is positive."""
    if pos <= 0:
        sys.stderr.write(f"Warning: Invalid {pos_name} for precursor {pri_id}: {pos}\n")
        return False
    return True


def _calculate_positive_strand_positions(excise_beg: int, mature_beg_rel: int, mature_end_rel: int,
                                        star_beg_rel: int, star_end_rel: int, mature_arm: str) -> Dict:
    """Calculate genomic positions for positive strand."""
    positions = {
        'mature_bg_abs': excise_beg + mature_beg_rel - 1,
        'mature_end_abs': excise_beg + mature_end_rel - 1,
        'star_bg_abs': excise_beg + star_beg_rel - 1,
        'star_end_abs': excise_beg + star_end_rel - 1,
        'pre_bg_abs': 0,
        'pre_end_abs': 0
    }
    
    if mature_arm == 'first':
        positions['mature_end_abs'] += 1
        positions['pre_bg_abs'] = positions['mature_bg_abs']
        positions['pre_end_abs'] = positions['star_end_abs']
    elif mature_arm == 'second':
        positions['mature_bg_abs'] -= 1
        positions['pre_bg_abs'] = positions['star_bg_abs'] - 1
        positions['pre_end_abs'] = positions['mature_end_abs']
    
    return positions


def _calculate_negative_strand_positions(excise_end: int, mature_beg_rel: int, mature_end_rel: int,
                                        star_beg_rel: int, star_end_rel: int, mature_arm: str) -> Dict:
    """Calculate genomic positions for negative strand."""
    positions = {
        'mature_bg_abs': excise_end - mature_end_rel,
        'mature_end_abs': excise_end - mature_beg_rel,
        'star_bg_abs': excise_end - star_end_rel,
        'star_end_abs': excise_end - star_beg_rel,
        'pre_bg_abs': 0,
        'pre_end_abs': 0
    }
    
    if mature_arm == 'first':
        positions['mature_bg_abs'] -= 1
        positions['pre_bg_abs'] = positions['star_bg_abs']
        positions['pre_end_abs'] = positions['mature_end_abs']
    elif mature_arm == 'second':
        positions['mature_end_abs'] += 1
        positions['pre_end_abs'] = positions['star_end_abs'] + 1
        positions['pre_bg_abs'] = positions['mature_bg_abs']
    
    return positions


def calculate_genomic_positions(predictions: List[Dict], precursors: Dict[str, Dict], 
                               chr_lengths: Dict[str, int]) -> List[Dict]:
    """Calculate genomic positions for predictions.
    
    Args:
        predictions: List of miRNA prediction dictionaries
        precursors: Dictionary of precursor information
        chr_lengths: Dictionary of chromosome lengths
        
    Returns:
        List of predictions with genomic positions added
    """
    results = []
    
    for pred in predictions:
        pri_id = pred.get('pri_id', '')
        if not pri_id:
            continue
        
        # precursor_info
        precursor_info = precursors.get(pri_id)
        if not precursor_info:
            found = False
            for key in precursors.keys():
                if key.startswith(pri_id):
                    precursor_info = precursors[key]
                    found = True
                    break
            
            if not found:
                for key in precursors.keys():
                    base_key = '_'.join(key.split('_')[:-1]) if '_' in key else key
                    if base_key == pri_id:
                        precursor_info = precursors[key]
                        found = True
                        break
            
            if not found:
                sys.stderr.write(f"Warning: Precursor {pri_id} not found in precursor file\n")
                continue
        
        # retieve relative positions
        try:
            mature_beg_rel = int(pred.get('mature_beg', 0))
            mature_end_rel = int(pred.get('mature_end', 0))
            star_beg_rel = int(pred.get('star_beg', 0))
            star_end_rel = int(pred.get('star_end', 0))
            mature_arm = pred.get('mature_arm', '')
        except (ValueError, TypeError) as e:
            sys.stderr.write(f"Warning: Invalid position values for precursor {pri_id}: {e}\n")
            continue
        
        strand = precursor_info['strand']
        chrom = precursor_info['chromosome']
        excise_beg = precursor_info['excise_beg']
        excise_end = precursor_info['excise_end']
        
        if chrom not in chr_lengths:
            sys.stderr.write(f"Warning: Chromosome {chrom} not found in chromosome length file\n")
            continue
        
        mature_bg_abs = 0
        mature_end_abs = 0
        star_bg_abs = 0
        star_end_abs = 0
        pre_bg_abs = 0
        pre_end_abs = 0
        
        if strand == '+':
            mature_bg_abs = excise_beg + mature_beg_rel - 1
            mature_end_abs = excise_beg + mature_end_rel - 1
            star_bg_abs = excise_beg + star_beg_rel - 1
            star_end_abs = excise_beg + star_end_rel - 1
            
            # Adjust precursor boundaries based on mature_arm
            if mature_arm == 'first':
                mature_bg_abs -= 1
                mature_end_abs -= 1
                pre_bg_abs = mature_bg_abs
                pre_end_abs = star_end_abs
            elif mature_arm == 'second':
                pre_bg_abs = star_bg_abs - 1
                pre_end_abs = mature_end_abs
            else:
                pre_bg_abs = excise_beg
                pre_end_abs = excise_end
            
            # Final adjustments
            if mature_arm == 'first':
                mature_end_abs += 1
            elif mature_arm == 'second':
                mature_bg_abs -= 1
                
        else:
            mature_bg_abs = excise_end - mature_end_rel
            mature_end_abs = excise_end - mature_beg_rel
            star_bg_abs = excise_end - star_end_rel
            star_end_abs = excise_end - star_beg_rel
            
            # Adjust precursor boundaries based on mature_arm
            if mature_arm == 'first':
                mature_bg_abs += 1
                mature_end_abs += 1
                pre_bg_abs = star_end_abs - 1
                pre_end_abs = mature_end_abs
            elif mature_arm == 'second':
                pre_bg_abs = mature_bg_abs
                pre_end_abs = star_end_abs
            else:
                pre_bg_abs = excise_beg
                pre_end_abs = excise_end
            
            # Final adjustments
            if mature_arm == 'first':
                mature_bg_abs -= 1
                pre_bg_abs = star_bg_abs
            elif mature_arm == 'second':
                mature_end_abs += 1
                pre_end_abs += 1
        
        # Verify coordinates are positive
        if (mature_bg_abs <= 0 or mature_end_abs <= 0 or 
            star_bg_abs <= 0 or star_end_abs <= 0 or
            pre_bg_abs <= 0 or pre_end_abs <= 0):
            sys.stderr.write(f"Warning: Invalid coordinates for precursor {pri_id}: negative or zero positions\n")
            continue
        
        if pre_bg_abs > pre_end_abs:
            pre_bg_abs, pre_end_abs = pre_end_abs, pre_bg_abs
        
        if mature_bg_abs > mature_end_abs:
            mature_bg_abs, mature_end_abs = mature_end_abs, mature_bg_abs
        
        if star_bg_abs > star_end_abs:
            star_bg_abs, star_end_abs = star_end_abs, star_bg_abs
        
        # Create result
        result = pred.copy()
        result.update({
            'chromosome': chrom,
            'strand': strand,
            'mature_bg_abs': mature_bg_abs,
            'mature_end_abs': mature_end_abs,
            'star_bg_abs': star_bg_abs,
            'star_end_abs': star_end_abs,
            'pre_bg_abs': pre_bg_abs,
            'pre_end_abs': pre_end_abs,
            'mature_beg_rel': mature_beg_rel,
            'mature_end_rel': mature_end_rel,
            'star_beg_rel': star_beg_rel,
            'star_end_rel': star_end_rel,
        })
        
        results.append(result)
    
    print(f"Calculated genomic positions for {len(results)} predictions")
    return results


def remove_redundant_predictions(predictions: List[Dict]) -> List[Dict]:
    """Remove redundant predictions with same chromosome and mature miRNA start.
    
    Args:
        predictions: List of miRNA prediction dictionaries
        
    Returns:
        List of non-redundant predictions
    """
    unique_predictions = {}
    
    for pred in predictions:
        chrom = pred.get('chromosome', '')
        mature_start = pred.get('mature_bg_abs', 0)
        key = (chrom, mature_start)
        
        if key not in unique_predictions:
            unique_predictions[key] = pred
    
    non_redundant = list(unique_predictions.values())
    print(f"Reduced from {len(predictions)} to {len(non_redundant)} non-redundant predictions")
    return non_redundant


def analyze_read_distribution(pred: Dict, total_reads: int) -> Dict:
    """Analyze read distribution for a prediction.
    
    Args:
        pred: miRNA prediction dictionary
        total_reads: Total number of reads
        
    Returns:
        Dictionary with read distribution statistics
    """
    mature_beg_rel = int(pred.get('mature_beg_rel', 0))
    mature_end_rel = int(pred.get('mature_end_rel', 0))
    star_beg_rel = int(pred.get('star_beg_rel', 0))
    star_end_rel = int(pred.get('star_end_rel', 0))
    
    mature_count = 0
    star_count = 0
    total_count = 0
    biggest_other = 0
    star_present = False
    
    if 'alignment_lines' in pred:
        for line in pred['alignment_lines']:
            parts = line.split('\t')
            if len(parts) < 6:
                continue
            
            read_id = parts[0]
            read_count_match = re.search(r'_x(\d+)$', read_id)
            if not read_count_match:
                continue
            
            read_count = int(read_count_match.group(1))
            total_count += read_count
            
            position_range = parts[5]
            pos_match = re.match(r'(\d+)\.\.(\d+)$', position_range)
            if not pos_match:
                continue
            
            pos_start = int(pos_match.group(1))
            pos_end = int(pos_match.group(2))
            
            if (star_beg_rel - 1 <= pos_start <= star_end_rel + 1 and 
                star_beg_rel - 1 <= pos_end <= star_end_rel + 1):
                star_count += read_count
                star_present = True
            elif (mature_beg_rel - 1 <= pos_start <= mature_end_rel + 1 and 
                  mature_beg_rel - 1 <= pos_end <= mature_end_rel + 1):
                mature_count += read_count
            else:
                if read_count > biggest_other:
                    biggest_other = read_count
    else:
        mature_query = pred.get('mature_query', '')
        if mature_query:
            match = re.search(r'_x(\d+)$', mature_query)
            if match:
                mature_count = int(match.group(1))
                total_count = mature_count
        
        star_seq = pred.get('star_seq', '')
        if star_seq and star_seq != 'None':
            star_present = True
    
    mature_rpm = (mature_count / total_reads * 1000000) if total_reads > 0 else 0
    
    return {
        'mature_count': mature_count,
        'star_count': star_count,
        'total_count': total_count,
        'biggest_other': biggest_other,
        'star_present': star_present,
        'mature_rpm': mature_rpm
    }


def apply_plant_criteria(predictions: List[Dict], total_reads: int) -> List[Dict]:
    """Apply plant-specific miRNA criteria to filter predictions.
    
    Args:
        predictions: List of miRNA prediction dictionaries
        total_reads: Total number of reads
        
    Returns:
        List of predictions meeting plant criteria
    """
    filtered_predictions = []
    
    print(f"Applying plant criteria to {len(predictions)} predictions...")
    print(f"Total reads: {total_reads}")
    
    for pred in predictions:
        mature_struct = pred.get('mature_struct', '')
        star_struct = pred.get('star_struct', '')
        
        if not mature_struct or not star_struct:
            continue
        
        mature_stem = mature_struct[:-2] if len(mature_struct) >= 2 else mature_struct
        star_stem = star_struct[:-2] if len(star_struct) >= 2 else star_struct
        
        mature_mismatch = mature_stem.count('.')
        star_mismatch = star_stem.count('.')
        
        if mature_mismatch > 5 or star_mismatch > 5 or abs(mature_mismatch - star_mismatch) > 3:
            continue
        
        read_stats = analyze_read_distribution(pred, total_reads)
        
        if read_stats['total_count'] == 0:
            continue
        
        mature_star_count = read_stats['mature_count'] + read_stats['star_count']
        if mature_star_count < 0.75 * read_stats['total_count']:
            continue
        
        mature_seq = pred.get('mature_seq', '')
        mature_len = len(mature_seq)
        
        if mature_len < 20 or mature_len > 24:
            continue
        
        if mature_len == 23 or mature_len == 24:
            if mature_mismatch == 0 or star_mismatch == 0:
                continue
            
            if read_stats['mature_rpm'] < 20:
                continue
            
            if not read_stats['star_present']:
                continue
        
        filtered_predictions.append(pred)
    
    print(f"Filtered from {len(predictions)} to {len(filtered_predictions)} predictions meeting plant criteria")
    return filtered_predictions


def write_predictions(predictions: List[Dict], output_file: str):
    """Write predictions to output file.
    
    Args:
        predictions: List of miRNA prediction dictionaries
        output_file: Path to output file
    """
    try:
        with open(output_file, 'w') as f:
            for pred in predictions:
                chrom = pred.get('chromosome', '')
                strand = pred.get('strand', '')
                mature_query = pred.get('mature_query', '')
                pri_id = pred.get('pri_id', '')
                
                mature_bg = pred.get('mature_bg_abs', 0)
                mature_end = pred.get('mature_end_abs', 0)
                mature_pos = f"{mature_bg}..{mature_end}"
                
                pre_bg = pred.get('pre_bg_abs', 0)
                pre_end = pred.get('pre_end_abs', 0)
                pre_pos = f"{pre_bg}..{pre_end}"
                
                mature_seq = pred.get('mature_seq', '')
                pre_seq = pred.get('pre_seq', '')
                
                pre_seq_length = len(pre_seq)
                pre_pos_length = pre_end - pre_bg
                
                if pre_seq_length != pre_pos_length:
                    sys.stderr.write(f"Warning: Sequence length ({pre_seq_length}) does not match position length ({pre_pos_length}) for {pri_id}\n")
                
                line = f"{chrom}\t{strand}\t{mature_query}\t{pri_id}\t{mature_pos}\t{pre_pos}\t{mature_seq}\t{pre_seq}\n"
                f.write(line)
        
        print(f"Wrote {len(predictions)} predictions to {output_file}")
        
    except Exception as e:
        sys.stderr.write(f"Error writing output file {output_file}: {e}\n")
        sys.exit(1)


def write_bed_file(predictions: List[Dict], output_file: str):
    """Write predictions in BED format.
    
    Args:
        predictions: List of miRNA prediction dictionaries
        output_file: Path to output BED file
    """
    try:
        with open(output_file, 'w') as f:
            for pred in predictions:
                chrom = pred.get('chromosome', '')
                strand = pred.get('strand', '')
                mature_query = pred.get('mature_query', '')
                pri_id = pred.get('pri_id', '')
                
                pre_bg = pred.get('pre_bg_abs', 0)
                pre_end = pred.get('pre_end_abs', 0)
                
                pre_bg_bed = pre_bg - 1
                pre_end_bed = pre_end
                
                mature_arm = pred.get('mature_arm', '')
                
                if strand == '+':
                    arm_type = '5' if mature_arm == 'first' else '3'
                else:
                    arm_type = '3' if mature_arm == 'first' else '5'
                
                name = f"{pri_id}_{arm_type}"
                line = f"{chrom}\t{pre_bg_bed}\t{pre_end_bed}\t{name}\t0\t{strand}\n"
                f.write(line)
        
        print(f"Wrote {len(predictions)} predictions to BED file {output_file}")
        
    except Exception as e:
        sys.stderr.write(f"Error writing BED file {output_file}: {e}\n")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Remove redundant miRNA predictions and apply plant-specific criteria',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script processes miRDeep prediction results to:
  1. Remove redundant predictions (same chromosome and mature start position)
  2. Apply plant-specific miRNA criteria

Input files:
  chr_length:     Chromosome length file (chromosome<TAB>length)
  precursor:      Precursor sequences in FASTA format
  miRDeep_result: miRDeep prediction results file
  total_reads:    File containing total reads count

Output files:
  nr_predictions:     Non-redundant predictions
  filtered_predictions: Predictions meeting plant criteria
  bed_file:           Predictions in BED format (optional)
        """
    )
    
    parser.add_argument('chr_length',
                       help='Chromosome length file')
    parser.add_argument('precursor',
                       help='Precursor sequences file (FASTA)')
    parser.add_argument('miRDeep_result',
                       help='miRDeep prediction results file')
    parser.add_argument('total_reads',
                       help='File containing total reads count')
    parser.add_argument('-n', '--nr-output', required=True,
                       help='Output file for non-redundant predictions')
    parser.add_argument('-f', '--filtered-output', required=True,
                       help='Output file for predictions meeting plant criteria')
    parser.add_argument('-b', '--bed-output',
                       help='Output file in BED format (optional)')
    parser.add_argument('--debug', action='store_true',
                       help='Enable debug output')
    parser.add_argument('--stats', action='store_true',
                       help='Show detailed statistics')
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("REMOVE REDUNDANT miRNA PREDICTIONS AND APPLY PLANT CRITERIA")
    print("=" * 70)
    print(f"Chromosome length file: {args.chr_length}")
    print(f"Precursor file:         {args.precursor}")
    print(f"miRDeep results:        {args.miRDeep_result}")
    print(f"Total reads file:       {args.total_reads}")
    print(f"Non-redundant output:   {args.nr_output}")
    print(f"Filtered output:        {args.filtered_output}")
    if args.bed_output:
        print(f"BED output:             {args.bed_output}")
    print("-" * 70)
    
    print("1. Parsing input files...")
    chr_lengths = parse_chromosome_lengths(args.chr_length)
    precursors = parse_precursors(args.precursor)
    predictions = parse_miRDeep_result(args.miRDeep_result)
    total_reads = get_total_reads(args.total_reads)
    
    print("2. Calculating genomic positions...")
    predictions_with_pos = calculate_genomic_positions(predictions, precursors, chr_lengths)
    
    if not predictions_with_pos:
        print("ERROR: No predictions with genomic positions calculated")
        sys.exit(1)
    
    print("3. Removing redundant predictions...")
    non_redundant = remove_redundant_predictions(predictions_with_pos)
    
    print("4. Writing non-redundant predictions...")
    write_predictions(non_redundant, args.nr_output)
    
    print("5. Applying plant-specific criteria...")
    filtered = apply_plant_criteria(non_redundant, total_reads)
    
    print("6. Writing filtered predictions...")
    write_predictions(filtered, args.filtered_output)
    
    if args.bed_output:
        print("7. Writing BED file...")
        write_bed_file(filtered, args.bed_output)
    
    if args.stats:
        print("\n" + "=" * 70)
        print("STATISTICS")
        print("=" * 70)
        print(f"Total predictions parsed:          {len(predictions)}")
        print(f"Predictions with positions:        {len(predictions_with_pos)}")
        print(f"Non-redundant predictions:         {len(non_redundant)}")
        print(f"Predictions meeting plant criteria: {len(filtered)}")
        
        if len(predictions) > 0:
            reduction = ((len(predictions)-len(non_redundant))/len(predictions)*100)
            print(f"Reduction rate (non-redundant):    {reduction:.1f}%")
        
        if len(non_redundant) > 0:
            filter_rate = ((len(non_redundant)-len(filtered))/len(non_redundant)*100)
            print(f"Filter rate (plant criteria):      {filter_rate:.1f}%")
        
        if len(predictions) > 0:
            success_rate = (len(filtered)/len(predictions)*100)
            print(f"Overall success rate:              {success_rate:.1f}%")
    
    print("\n" + "=" * 70)
    print("PROCESSING COMPLETE")
    print("=" * 70)
    print(f"Non-redundant predictions written to: {args.nr_output}")
    print(f"Filtered predictions written to:      {args.filtered_output}")
    if args.bed_output:
        print(f"BED file written to:                {args.bed_output}")
    print("=" * 70)


if __name__ == "__main__":
    main()
