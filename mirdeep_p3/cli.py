#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""Console-script entry point.

Locates the MirDeep-P3 checkout and runs its entry script in-process, so the
banner, argument parsing and PROJECT_ROOT behave exactly as when the script is
executed directly.
"""

import os
import runpy
import sys
from pathlib import Path

SCRIPT_NAME = "mirdeep-p3"


def find_project_root() -> Path:
    """Return the checkout that holds the mirdeep-p3 entry script.

    MIRDEEP_P3_HOME wins; otherwise the editable-install layout (this package
    next to the entry script) and the current directory are tried.
    """
    env = os.environ.get("MIRDEEP_P3_HOME")
    if env:
        root = Path(env).expanduser().resolve()
        if (root / SCRIPT_NAME).is_file():
            return root
        raise SystemExit(f"mirdeep-p3: MIRDEEP_P3_HOME={env} has no {SCRIPT_NAME} script")

    here = Path(__file__).resolve().parent
    for root in (here.parent, here.parent.parent, Path.cwd(), Path.cwd().parent):
        if (root / SCRIPT_NAME).is_file() and (root / "src" / "version.py").is_file():
            return root.resolve()
    raise SystemExit(
        "mirdeep-p3: cannot locate the MirDeep-P3 checkout; "
        "set MIRDEEP_P3_HOME to the directory that contains the mirdeep-p3 script")


def main(argv=None) -> None:
    root = find_project_root()
    script = root / SCRIPT_NAME
    if argv is None:
        sys.argv = [str(script)] + sys.argv[1:]
    else:
        sys.argv = [str(script)] + list(argv)
    runpy.run_path(str(script), run_name="__main__")


if __name__ == "__main__":
    main()
