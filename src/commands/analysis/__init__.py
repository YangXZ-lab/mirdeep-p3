"""
MirDeep-P3 analysis sub-package.
Contains subcommands: TFBS, Target_finder, Differential_expression,
                     Functional_analysis, Stat, Onestep.
"""

import sys
import argparse

from . import (tfbs, target_finder, differential_expression,
               functional_analysis, stat, onestep)

def main(argv, project_root=None):
    """Entry point for analysis subcommands."""
    parser = argparse.ArgumentParser(
        prog="mirdeep-p3.py analysis",
        description="Downstream analysis modules for MirDeep-P3.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(
        title="analysis subcommands",
        dest="analysis_command",
        required=True,
    )

    # TFBS
    parser_tfbs = subparsers.add_parser(
        "TFBS",
        help="Transcription factor binding site analysis.",
    )
    tfbs.add_arguments(parser_tfbs)

    # Target_finder
    parser_target = subparsers.add_parser(
        "Target_finder",
        help="miRNA target prediction.",
    )
    target_finder.add_arguments(parser_target)

    # Differential_expression
    parser_de = subparsers.add_parser(
        "Differential_expression",
        aliases=["DE"],
        help="Differential expression analysis.",
        add_help=False,
    )
    differential_expression.add_arguments(parser_de)

    # Functional_analysis
    parser_func = subparsers.add_parser(
        "Functional_analysis",
        aliases=["FA"],
        help="Functional enrichment analysis.",
    )
    functional_analysis.add_arguments(parser_func)

     # Stat
    parser_stat = subparsers.add_parser(
        "Stat",
        help="Basic statistics and structure plots from annotation results.",
        add_help=False,
    )
    stat.add_arguments(parser_stat)
    
    # Onestep
    parser_onestep = subparsers.add_parser(
        "Onestep",
        help="Run full analysis pipeline from basic-info to functional enrichment.",
        add_help=False,
    )
    onestep.add_arguments(parser_onestep)
    
    args = parser.parse_args(argv)
    args.project_root = project_root   # 传递给子命令

    # A subparser reached through an alias reports the alias itself in
    # analysis_command ("DE", "FA"), so normalise the name before dispatching;
    # otherwise the alias falls through to the help branch and exits 1.
    command = {"DE": "Differential_expression",
               "FA": "Functional_analysis"}.get(args.analysis_command,
                                                args.analysis_command)

    if command == "TFBS":
        tfbs.run(args)
    elif command == "Target_finder":
        target_finder.run(args)
    elif command == "Differential_expression":
        differential_expression.run(args)
    elif command == "Functional_analysis":
        functional_analysis.run(args)
    elif command == "Stat":
        stat.run(args)
    elif command == "Onestep":
        onestep.run(args)
    else:
        parser.print_help()
        sys.exit(1)
