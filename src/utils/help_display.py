#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""Combined help for runs that select more than one step (--step 1,2 / 1,3 ...).

The text is generated from the steps' own parsers instead of being written out
by hand, so it cannot drift away from the options that actually exist: a new
argument shows up here automatically, and no combination falls back to a
placeholder any more (it used to for everything except --step 1,2).
"""

import argparse
import sys
import textwrap
from typing import Dict, List, Sequence

# Canonical order; any step not listed here is appended in user order.
STEP_ORDER = ("identification", "annotation", "analysis")

_WIDTH = 100
_INDENT = "  "
_GAP = "    "
_MAX_LEFT = 42


def _is_help_action(action) -> bool:
    return action.dest == "help" or set(action.option_strings) == {"-h", "--help"}


def _invocation(action) -> str:
    """Render an option the way the usage line does, e.g. '-i INPUT, --input INPUT'."""
    parts = []
    for option in action.option_strings:
        if action.nargs == 0:
            parts.append(option)
            continue
        metavar = action.metavar or action.dest.upper()
        if isinstance(metavar, tuple):
            metavar = metavar[0]
        parts.append(f"{option} {metavar}")
    return ", ".join(parts)


def _print_rows(out, entries: Sequence) -> None:
    """Print invocation + wrapped help, aligned in two columns."""
    invocations = [_invocation(entry["action"]) for entry in entries]
    left = min(max(len(text) for text in invocations) + len(_GAP), _MAX_LEFT)
    body = max(20, _WIDTH - len(_INDENT) - left)
    for entry, invocation in zip(entries, invocations):
        help_text = entry["helps"][0]
        if len(set(entry["helps"])) > 1:
            help_text += "   (help text differs per step)"
        if not help_text:
            out.write(f"{_INDENT}{invocation}\n")
            continue
        lines = textwrap.wrap(help_text, body)
        out.write(f"{_INDENT}{invocation.ljust(left)}{lines[0]}\n")
        for line in lines[1:]:
            out.write(f"{_INDENT}{' ' * left}{line}\n")


def _subcommand_help(action) -> List:
    """Return [(name, help)] for a subparsers action, skipping alias entries."""
    helps = {choice.dest: " ".join((choice.help or "").split())
             for choice in getattr(action, "_choices_actions", [])}
    seen = set()
    entries = []
    for name, sub in action.choices.items():
        if id(sub) in seen:      # an alias of an already listed subcommand
            continue
        seen.add(id(sub))
        entries.append((name, helps.get(name, "")))
    return entries


def _collect(steps: Sequence[str], parsers: Dict[str, argparse.ArgumentParser]):
    """Group every option by the set of steps that accepts it."""
    groups: Dict[tuple, List] = {}
    subcommands: Dict[str, List] = {}
    present: List[str] = []

    for step in steps:
        parser = parsers.get(step)
        if parser is None:
            continue
        present.append(step)
        for action in parser._actions:
            if _is_help_action(action):
                continue
            if isinstance(action, argparse._SubParsersAction):
                subcommands[step] = _subcommand_help(action)
                continue
            if not action.option_strings:
                continue
            key = tuple(action.option_strings)
            entry = groups.setdefault(key, {"action": action, "steps": [],
                                            "helps": []})
            entry["steps"].append(step)
            entry["helps"].append(" ".join((action.help or "").split()))

    owners: Dict[tuple, List] = {}
    for entry in groups.values():
        accepted = tuple(s for s in present if s in entry["steps"])
        owners.setdefault(accepted, []).append(entry)
    return present, owners, subcommands


def print_combined_step_help(steps: Sequence[str],
                             parsers: Dict[str, argparse.ArgumentParser]) -> None:
    """Print generated combined help for the given steps.

    Args:
        steps: canonical step names, in the order they should be shown.
        parsers: step name -> the parser that step uses.
    """
    out = sys.stdout
    present, owners, subcommands = _collect(steps, parsers)
    if not present:
        out.write("No recognised steps to describe.\n")
        return

    width = 62
    out.write("\n" + "\u2554" + "\u2550" * width + "\u2557\n")
    out.write("\u2551" + "  MirDeep-P3  combined help".ljust(width) + "\u2551\n")
    out.write("\u255a" + "\u2550" * width + "\u255d\n")
    out.write(f"Steps: {', '.join(present)}\n")

    if "identification" in present and "annotation" in present:
        out.write("\nCombined run layout:\n")
        out.write("  OUT/identification/   identification results (prefix subfolders + .pipe)\n")
        out.write("  OUT/annotation/       annotation results\n")
        out.write("  annotation -i is filled from OUT/identification/mirdp3-identification-*.pipe\n")
        out.write("  --prefix_miRNA and --species are required for annotation;\n")
        out.write("  -g / -d / -r / -t / -p fall back to the identification values.\n")

    ordered = sorted(
        owners.items(),
        key=lambda item: (len(item[0]) != len(present),
                          [present.index(s) for s in item[0]],
                          -len(item[1])))

    for accepted, actions in ordered:
        if len(accepted) == len(present):
            title = "Accepted by every selected step"
        else:
            title = "Accepted by: " + ", ".join(accepted)
        out.write(f"\n{title}\n" + "-" * len(title) + "\n")
        _print_rows(out, actions)

    for step, entries in subcommands.items():
        if not entries:
            continue
        out.write(f"\n{step} subcommands\n" + "-" * (len(step) + 12) + "\n")
        left = min(max(len(name) for name, _ in entries) + len(_GAP), _MAX_LEFT)
        for name, help_text in entries:
            lines = textwrap.wrap(help_text, max(20, _WIDTH - len(_INDENT) - left)) or [""]
            out.write(f"{_INDENT}{name.ljust(left)}{lines[0]}\n")
            for line in lines[1:]:
                out.write(f"{_INDENT}{' ' * left}{line}\n")

    out.write("\nRun 'mirdeep-p3 <step> -h' for the complete help of a single step.\n")
