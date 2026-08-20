#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
Check all required Python dependencies for MirDeep-P3.
Add third-party packages to REQUIRED_PACKAGES as needed.
"""

import sys
import importlib

REQUIRED_PACKAGES = [
    "sys",
    "os",
    "gzip",
    "argparse",
    "typing",
    "collections",
    "re",
    "pathlib",
    "importlib",
    "subprocess",
    "tempfile",
    "traceback"
]

def check_python_deps():
    missing = []
    print("Python package dependency check:")
    if not REQUIRED_PACKAGES:
        print("  (no third-party dependencies required)")
    for pkg in REQUIRED_PACKAGES:
        try:
            importlib.import_module(pkg)
            print(f"  {pkg} ✔")
        except ImportError:
            print(f"  {pkg} ✘")
            missing.append(pkg)
    if missing:
        print("\nERROR: The following Python packages are missing:")
        for m in missing:
            print(f"  - {m}")
        sys.exit(1)
    else:
        print("All Python dependencies satisfied.")

if __name__ == "__main__":
    check_python_deps()
