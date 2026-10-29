#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
Data-index release bookkeeping (O-10).

The data index (Rfam / mature / PmiREN Blast databases) and the software are
published on separate schedules, and the software is updated far more often than
the index.  Comparing the index with the software version therefore says
nothing useful.  What is worth knowing is which index release is installed, and
whether it is the one data-index.env currently names -- so this module compares
"installed index tag" with "expected index tag", both on the data axis.

The tag is written to <data>/index/.index-version when the index is unpacked
from data-index.tar.gz (test.sh does it, and packagers should).  When the file
is absent the state is "unknown": that is not an error, the code just cannot
tell which release is on disk.

Nothing here changes a result or aborts a run; the information goes to the
console and into the run manifest.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

ENV_FILE = "data-index.env"
VERSION_FILE = ".index-version"
ENV_KEY = "DATA_INDEX_RELEASE"
INDEX_SUBDIR = Path("data") / "index"


def expected_release(project_root) -> Optional[str]:
    """Return DATA_INDEX_RELEASE from data-index.env, or None if unavailable."""
    try:
        text = (Path(project_root) / ENV_FILE).read_text(encoding="utf-8")
    except OSError:
        return None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(ENV_KEY + "="):
            value = line.split("=", 1)[1].strip().strip("'\"")
            return value or None
    return None


def installed_release(index_dir) -> Optional[str]:
    """Return the tag recorded in <index_dir>/.index-version, or None."""
    try:
        text = (Path(index_dir) / VERSION_FILE).read_text(encoding="utf-8").strip()
    except OSError:
        return None
    return text or None


def write_installed(index_dir, release: str) -> bool:
    """Record `release` as the installed tag. Returns False if it cannot be written."""
    target = Path(index_dir) / VERSION_FILE
    try:
        target.write_text(str(release).strip() + "\n", encoding="utf-8")
    except OSError:
        return False
    return True


def index_status(project_root, index_dir=None) -> Dict[str, Optional[str]]:
    """Describe the installed index relative to data-index.env.

    state is one of:
        missing    -- no index directory at all
        unknown    -- index present but no .index-version marker
        unrecorded -- index and marker present, data-index.env has no tag
        match      -- installed tag equals DATA_INDEX_RELEASE
        mismatch   -- installed tag differs from DATA_INDEX_RELEASE
    """
    root = Path(project_root)
    index = Path(index_dir) if index_dir else root / INDEX_SUBDIR
    expected = expected_release(root)
    installed = installed_release(index)

    if not index.is_dir():
        state = "missing"
    elif installed is None:
        state = "unknown"
    elif expected is None:
        state = "unrecorded"
    elif installed == expected:
        state = "match"
    else:
        state = "mismatch"

    return {
        "state": state,
        "installed": installed,
        "expected": expected,
        "path": str(index),
    }


def describe_index(status: Dict[str, Optional[str]]) -> str:
    """One-line, human-readable form of an index_status() result."""
    state = status.get("state")
    installed = status.get("installed")
    expected = status.get("expected")
    path = status.get("path", "")

    if state == "missing":
        return f"not found at {path} (identification and annotation need it)"
    if state == "match":
        return f"{installed} (matches {ENV_FILE})"
    if state == "mismatch":
        return (f"{installed} differs from {ENV_FILE} ({expected}); the run is fine, "
                f"but results may differ from one made with the published index")
    if state == "unknown":
        return (f"release unknown at {path} (no {VERSION_FILE}; "
                f"{ENV_FILE} expects {expected})")
    if state == "unrecorded":
        return f"{installed} (no tag in {ENV_FILE})"
    return f"state {state!r} at {path}"
