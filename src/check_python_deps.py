#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
Python dependency check for MirDeep-P3 (I-04).

MirDeep-P3 itself uses only the standard library.  The previous version of this
file listed a dozen stdlib modules under REQUIRED_PACKAGES, so its "All Python
dependencies satisfied" line proved nothing.  It now reports what is actually
required, and treats the one real third-party package as optional because it is
only imported by one engine:

    numpy -- src/MiRNATarget/psrna_scan.py, i.e. Target_finder --engine psrna
             (the default engine).  identification and annotation never touch it.

Exit status is 1 only when a required package is missing.
"""

import importlib
import sys

# Nothing outside the standard library is needed to run the pipeline.
REQUIRED_PACKAGES = []

# Optional: key -> (import name, why it is needed)
OPTIONAL_PACKAGES = {
    "numpy": ("numpy", "Target_finder --engine psrna (src/MiRNATarget/psrna_scan.py)"),
}


def _status(name):
    try:
        module = importlib.import_module(name)
    except ImportError:
        return None
    return getattr(module, "__version__", "")


def check_python_deps():
    print("Python dependency check:")
    missing = []
    if not REQUIRED_PACKAGES:
        print("  (no third-party dependencies required to run the pipeline)")
    for pkg in REQUIRED_PACKAGES:
        version = _status(pkg)
        if version is None:
            print("  %s \u2718" % pkg)
            missing.append(pkg)
        else:
            print("  %s \u2714 %s" % (pkg, version))

    for label, (module_name, why) in OPTIONAL_PACKAGES.items():
        version = _status(module_name)
        if version is None:
            print("  %s \u2718 optional -- %s will not work" % (label, why))
        else:
            print("  %s \u2714 %s (optional: %s)" % (label, version, why))

    print("  python %d.%d.%d" % sys.version_info[:3])
    if missing:
        print("\nERROR: The following Python packages are missing:")
        for name in missing:
            print("  - %s" % name)
        sys.exit(1)
    print("All required Python dependencies satisfied.")


if __name__ == "__main__":
    check_python_deps()
