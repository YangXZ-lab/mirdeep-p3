#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
Annotate prediction BED file with conservation status based on cluster file.

"""

import sys
import re
import argparse
from typing import Dict


def load_cluster_file(cluster_file: str, strip_strand: bool = True) -> Dict[str, str]:
    """
    Load cluster file and build a mapping from name to family.

    Cluster file format (tab-separated):
        col1: name with strand suffix e.g. "N710-Chr01_76-...-4(+)"
        col2: family, e.g. "MIR156" or "MIRN1"
        (additional columns, if any, are ignored)

    Args:
        cluster_file: Path to cluster file.
        strip_strand: If True, remove trailing "(+)" or "(-)" from col1.

    Returns:
        Dictionary mapping name -> family.
    """
    mapping = {}
    try:
        with open(cluster_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    print(f"Warning: Line {line_num} in {cluster_file} "
                          f"has fewer than 2 columns, skipping: {line}",
                          file=sys.stderr)
                    continue
                name_with_strand = parts[0].strip()
                family = parts[1].strip()

                if strip_strand:
                    # Remove trailing strand annotation like (+), (-), (.), etc.
                    name = re.sub(r'\([+-]\)$', '', name_with_strand).strip()
                else:
                    name = name_with_strand

                if name in mapping:
                    print(f"Warning: Duplicate name '{name}' in cluster file "
                          f"(line {line_num}), overwriting previous mapping.",
                          file=sys.stderr)
                mapping[name] = family
    except Exception as e:
        print(f"Error: Cannot read cluster file {cluster_file}: {e}", file=sys.stderr)
        sys.exit(1)

    if not mapping:
        print(f"Warning: No valid entries found in cluster file {cluster_file}",
              file=sys.stderr)
    return mapping


def annotate_bed_file(bed_file: str, cluster_mapping: Dict[str, str],
                      output_file: str) -> None:
    """
    Annotate each line of the BED file with conservation status.

    For each BED entry:
      - If name matches a cluster entry, and the family starts with 'MIR' but
        not 'MIRN', annotate as 'conserved=<family>'.
      - Otherwise (including MIRN families or no match), annotate as
        'non_conserved'.

    Args:
        bed_file: Input BED file (8 columns, tab-separated).
        cluster_mapping: Dictionary mapping name -> family.
        output_file: Output file path (original 8 columns + annotation column).
    """
    processed = 0
    conserved = 0
    non_conserved = 0
    skipped = 0

    try:
        with open(bed_file, 'r') as infh, open(output_file, 'w') as outfh:
            for line_num, line in enumerate(infh, 1):
                line = line.rstrip('\n')
                if not line:
                    continue

                fields = line.split('\t')
                if len(fields) < 8:
                    print(f"Warning: Line {line_num} in {bed_file} "
                          f"has fewer than 8 columns, skipping: {line}",
                          file=sys.stderr)
                    skipped += 1
                    continue

                # The 4th column (index 3) is the name
                name = fields[3].strip()

                family = cluster_mapping.get(name)
                if family is not None and family.startswith('MIR') and not family.startswith('MIRN'):
                    annotation = f"conserved={family}"
                    conserved += 1
                else:
                    annotation = "non_conserved"
                    non_conserved += 1

                outfh.write(f"{line}\t{annotation}\n")
                processed += 1
    except Exception as e:
        print(f"Error: Cannot process file {bed_file}: {e}", file=sys.stderr)
        sys.exit(1)

    # Summary
    print(f"Annotation completed.", file=sys.stderr)
    print(f"  Total processed lines: {processed}", file=sys.stderr)
    print(f"  Conserved (known MIR): {conserved}", file=sys.stderr)
    print(f"  Non-conserved (MIRN or unmatched): {non_conserved}", file=sys.stderr)
    if skipped:
        print(f"  Skipped lines: {skipped}", file=sys.stderr)
    print(f"  Output written to: {output_file}", file=sys.stderr)


def main():
    """Command line interface."""
    parser = argparse.ArgumentParser(
        description='Annotate prediction BED file with conservation status '
                    'using a cluster file.',
        epilog='The cluster file should be tab-separated with at least two columns: '
               'name_with_strand, family. '
               'By default, strand suffix (e.g., "(+)") is stripped from the name '
               'before matching.',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '-i', '--input-bed',
        required=True,
        help='Input BED file (8 columns, tab-separated).'
    )

    parser.add_argument(
        '-c', '--cluster',
        required=True,
        help='Cluster file (tab-separated: name_with_strand, family).'
    )

    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output file (original 8 columns + annotation column).'
    )

    parser.add_argument(
        '--no-strip-strand',
        action='store_true',
        help='Do not strip strand suffix (e.g., "(+)") from cluster file names.'
    )

    args = parser.parse_args()

    strip_strand = not args.no_strip_strand
    cluster_mapping = load_cluster_file(args.cluster, strip_strand=strip_strand)

    if not cluster_mapping:
        print("Warning: Cluster mapping is empty. All entries will be annotated as 'non_conserved'.",
              file=sys.stderr)

    annotate_bed_file(args.input_bed, cluster_mapping, args.output)


if __name__ == "__main__":
    main()