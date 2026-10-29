# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""Launcher package for MirDeep-P3.

The pipeline itself lives in the checkout (the ``mirdeep-p3`` entry script,
``src/``, ``scripts/`` and ``data/``).  This package only provides the
``mirdeep-p3`` console script, so an editable install does not have to copy
(and then silently desynchronise) the tree.  Set ``MIRDEEP_P3_HOME`` when the
checkout is not the one this package was installed from.
"""

import os
from pathlib import Path


def _read_version() -> str:
    """Read __version__ from the checkout's src/version.py, the single source."""
    here = Path(__file__).resolve().parent
    candidates = []
    env = os.environ.get("MIRDEEP_P3_HOME")
    if env:
        candidates.append(Path(env).expanduser())
    candidates.extend([here.parent, here.parent.parent])
    for root in candidates:
        version_py = root / "src" / "version.py"
        try:
            text = version_py.read_text(encoding="utf-8")
        except OSError:
            continue
        for line in text.splitlines():
            if line.startswith("__version__"):
                return line.split("=", 1)[1].strip().strip("\"'")
    return "0.0.0"


__version__ = _read_version()
