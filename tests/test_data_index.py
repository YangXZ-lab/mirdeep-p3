#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
Unit tests for src/utils/data_index.py (O-10).

Standard library only: run with `python3 tests/test_data_index.py` or pytest.
The point of the module is that the index tag is compared with data-index.env,
never with the software version, so the tests pin every state.
"""

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from utils.data_index import (  # noqa: E402
    describe_index,
    expected_release,
    index_status,
    installed_release,
    write_installed,
)

TAG = "mirdeep-p3-v3.1.4c-full"


def _project(tmp: Path, env_line="DATA_INDEX_RELEASE=" + TAG + "\n", make_index=True):
    root = tmp
    (root / "data-index.env").write_text("# comment\n" + env_line, encoding="utf-8")
    if make_index:
        (root / "data" / "index").mkdir(parents=True)
    return root


def test_expected_release_reads_env():
    with tempfile.TemporaryDirectory() as td:
        root = _project(Path(td))
        assert expected_release(root) == TAG
        # Quotes and trailing spaces are tolerated.
        (root / "data-index.env").write_text('  DATA_INDEX_RELEASE="%s"  \n' % TAG,
                                             encoding="utf-8")
        assert expected_release(root) == TAG
        # A file without the key yields None.
        (root / "data-index.env").write_text("# nothing here\n", encoding="utf-8")
        assert expected_release(root) is None
        # A missing env file yields None.
        (root / "data-index.env").unlink()
        assert expected_release(root) is None


def test_states():
    with tempfile.TemporaryDirectory() as td:
        root = _project(Path(td))
        index = root / "data" / "index"

        # Index present, no marker -> unknown.
        status = index_status(root)
        assert status["state"] == "unknown", status
        assert status["installed"] is None
        assert status["expected"] == TAG
        assert "release unknown" in describe_index(status)

        # Marker matching the env -> match.
        assert write_installed(index, TAG)
        status = index_status(root)
        assert status["state"] == "match", status
        assert "matches data-index.env" in describe_index(status)

        # Marker differing from the env -> mismatch, and the message says results
        # may differ rather than claiming the run will fail.
        write_installed(index, "data-v2029")
        status = index_status(root)
        assert status["state"] == "mismatch", status
        message = describe_index(status)
        assert "differs from data-index.env" in message
        assert "may differ" in message

        # env without a tag -> unrecorded.
        (root / "data-index.env").write_text("# empty\n", encoding="utf-8")
        status = index_status(root)
        assert status["state"] == "unrecorded", status
        assert "no tag in data-index.env" in describe_index(status)


def test_missing_index_directory():
    with tempfile.TemporaryDirectory() as td:
        root = _project(Path(td), make_index=False)
        status = index_status(root)
        assert status["state"] == "missing", status
        assert "not found" in describe_index(status)


def test_explicit_index_dir_and_blank_marker():
    with tempfile.TemporaryDirectory() as td:
        base = Path(td)
        root = _project(base, make_index=False)
        other = base / "elsewhere"
        other.mkdir()
        # An explicitly passed directory wins over <project>/data/index.
        assert index_status(root, other)["state"] == "unknown"
        # A blank marker counts as unknown, not as an empty tag.
        (other / ".index-version").write_text("   \n", encoding="utf-8")
        assert installed_release(other) is None
        assert write_installed(other, TAG)
        assert installed_release(other) == TAG
        assert index_status(root, other)["state"] == "match"


def test_tag_is_data_axis_only():
    # The state must not depend on src/version.py at all.
    with tempfile.TemporaryDirectory() as td:
        root = _project(Path(td))
        (root / "src").mkdir()
        (root / "src" / "version.py").write_text('__version__ = "9.9.9"\n', encoding="utf-8")
        write_installed(root / "data" / "index", TAG)
        status = index_status(root)
        assert status["state"] == "match", status


def main() -> int:
    tests = [value for name, value in sorted(globals().items())
             if name.startswith("test_") and callable(value)]
    failures = 0
    for test in tests:
        try:
            test()
        except AssertionError as exc:
            failures += 1
            print(f"FAIL {test.__name__}: {exc}")
        else:
            print(f"ok   {test.__name__}")
    print(f"\n{len(tests) - failures}/{len(tests)} passed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
