#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
External-tool resolution for the analysis subcommands (I-03 / O-04).

Every analysis step used to spell the tool name straight into the shell command
("Rscript", "fimo", "ssearch36", ...).  A missing tool then surfaced as a bare
"command not found" from somewhere inside a pipeline, and a non-default build
could not be pinned.  Every subcommand now has the same shape:

    add_tool_arguments(parser)                       # --rscript/--fimo/...
    tools = resolve_tools(args, required=["fimo"])   # check once, up front

resolve_tools() returns a plain dict {key: path}.  An explicit command-line
value always wins; otherwise the tool is looked up in PATH.  Required tools are
exercised once (--version / -version / no argument) before the pipeline starts,
and a missing one aborts with the exact flag to use instead of failing later
inside a half-finished run.

Backward compatibility:
  * every flag is optional -- without them the historical PATH lookup is used;
  * the printed block keeps the format of check_external_tools(), so existing
    logs read the same.
"""

import argparse
import os
import shutil
import subprocess
import sys
from typing import Dict, Iterable, List, Sequence, Tuple

# key -> (default executable, CLI flag, what MirDeep-P3 uses it for)
TOOL_SPECS: Dict[str, Tuple[str, str, str]] = {
    "rscript": ("Rscript", "--rscript",
                "R scripts under scripts/ (statistics, PCA, DE, enrichment)"),
    "fimo": ("fimo", "--fimo", "MEME-suite FIMO motif scanning (TFBS)"),
    "bedtools": ("bedtools", "--bedtools",
                 "getfasta for upstream/promoter extraction (TFBS)"),
    "ssearch36": ("ssearch36", "--ssearch36",
                  "FASTA-suite ssearch36 (Target_finder --engine ssearch)"),
    "seqkit": ("seqkit", "--seqkit",
               "reverse complement and length bookkeeping (Target_finder --engine ssearch)"),
    "emapper": ("emapper.py", "--emapper",
                "eggNOG-mapper protein annotation (Functional_analysis -p)"),
    # 'python' has no PATH entry of its own: it defaults to the interpreter that
    # is running MirDeep-P3 (see resolve_tools).
    "python": ("python3", "--python",
               "helper scripts under scripts/ and src/ (default: this interpreter)"),
}


def add_tool_arguments(parser: argparse.ArgumentParser) -> None:
    """Add one --<tool> PATH override per entry in TOOL_SPECS."""
    group = parser.add_argument_group(
        "external tool paths",
        "Optional: pin a specific build instead of using whatever is in PATH.")
    for key, (default, flag, desc) in TOOL_SPECS.items():
        group.add_argument(flag, dest="tool_%s" % key, metavar="PATH",
                           help="%s (default: %s from PATH)" % (desc, default))


def _executable(path: str) -> bool:
    return bool(path) and os.path.isfile(path) and os.access(path, os.X_OK)


def _runs(command: str) -> bool:
    """Exercise the tool once, mirroring check_external_tools()."""
    for extra in ("--version", "-version", None):
        try:
            proc = subprocess.run([command] + ([extra] if extra else []),
                                  stdout=subprocess.DEVNULL,
                                  stderr=subprocess.DEVNULL, timeout=20)
        except (subprocess.TimeoutExpired, FileNotFoundError, PermissionError, OSError):
            return False
        if proc.returncode == 0:
            return True
    return False


def resolve_tools(args,
                  required: Sequence[str] = (),
                  optional: Sequence[str] = (),
                  quiet: bool = False,
                  hint: str = None) -> Dict[str, str]:
    """Resolve external tool paths and abort if a required one is unusable.

    args      : the parsed namespace (may carry tool_<key> overrides)
    required  : keys that must be present and runnable
    optional  : keys that are only needed on some code paths (e.g. --picture);
                a missing optional tool is simply absent from the result
    quiet     : suppress the status block (still aborts on a missing tool)
    hint      : extra line printed with the failure message
    """
    resolved: Dict[str, str] = {}
    missing: List[str] = []
    if not quiet:
        print("External dependency check:")

    for key in list(required) + [k for k in optional if k not in required]:
        if key == "python":
            continue  # handled below, always available
        default, flag, _ = TOOL_SPECS[key]
        is_required = key in required
        override = getattr(args, "tool_%s" % key, None)
        command = ""
        note = ""

        if override:
            candidate = os.path.expanduser(str(override))
            if _executable(candidate) and _runs(candidate):
                command = candidate
            else:
                note = "%s %s is not a working executable" % (flag, candidate)
        else:
            found = shutil.which(default)
            if found and _runs(found):
                command = found
            elif found:
                note = "found at %s but it did not run" % found
            else:
                note = "not found in PATH"

        if command:
            resolved[key] = command
            if not quiet:
                print("  %s \u2714 (%s)" % (default, command))
        else:
            if is_required or override:
                if not quiet:
                    print("  %s \u2718 (%s)" % (default, note))
            if is_required:
                missing.append(key)

    # The interpreter running MirDeep-P3 is always valid; only mention it when
    # the user pinned a different one.
    py_override = getattr(args, "tool_python", None)
    if py_override:
        candidate = os.path.expanduser(str(py_override))
        if not (_executable(candidate) and _runs(candidate)):
            sys.exit("Error: --python %s is not a working executable" % candidate)
        resolved["python"] = candidate
        if not quiet:
            print("  python \u2714 (%s)" % candidate)
    else:
        resolved["python"] = sys.executable or "python3"

    if missing:
        names = ", ".join(TOOL_SPECS[k][0] for k in missing)
        flags = "\n".join("  %s /path/to/%s" % (TOOL_SPECS[k][1], TOOL_SPECS[k][0])
                          for k in missing)
        message = ["Error: missing required external tool(s): %s" % names,
                   "Install them, activate the environment that provides them, or point to",
                   "them explicitly:", flags]
        if hint:
            message.append(hint)
        sys.exit("\n".join(message))

    return resolved


def forward_tool_args(args) -> List[str]:
    """Return [flag, path, ...] for every override set on args (used by Onestep)."""
    forwarded: List[str] = []
    for key, (default, flag, _) in TOOL_SPECS.items():
        value = getattr(args, "tool_%s" % key, None)
        if value:
            forwarded += [flag, str(value)]
    return forwarded
