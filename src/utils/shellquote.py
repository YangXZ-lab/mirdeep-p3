#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
Safe quoting of paths interpolated into shell command strings.

Several pipeline steps build a command line as a plain string and run it with
shell=True, because they rely on redirections and pipes.  Every path going into
such a string must be quoted: otherwise a name containing a space, a parenthesis,
a quote or any shell metacharacter is re-parsed by the shell, the command splits
into extra arguments, and the tool fails with a confusing non-zero status
(for example bowtie exiting 2).

Use it as::

    from utils.shellquote import shq
    cmd = f"{shq(tool)} -f {shq(input_fa)} > {shq(out)} 2>> {shq(log)}"

Only paths/names/values should be passed through shq() -- never shell syntax such
as ">", "2>&1", "|", and never a literal option string that is meant to expand into
several words.
"""

import shlex

__all__ = ["shq"]


def shq(value) -> str:
    """Return *value* quoted so that the shell treats it as a single word."""
    return shlex.quote(str(value))
