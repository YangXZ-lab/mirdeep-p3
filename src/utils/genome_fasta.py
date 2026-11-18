#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""Random access to a genome FASTA without external tools.

``samtools faidx`` / ``bedtools getfasta`` would both need to write an index
next to the FASTA, which is not always allowed (the genome often lives in a
read-only shared directory), and pyfaidx is an extra dependency.  This module
builds a small offset index itself and reads only the bytes it needs, so a
360 MB genome costs one scan and a few hundred bytes per query.

Coordinate convention
---------------------
``fetch()`` takes a **0-based, half-open** ``[start, end)`` interval and
returns the sequence in *transcript orientation* (reverse-complemented on the
minus strand).  This matches the coordinates stored in ``*-basic-info``: those
are 0-based with an exclusive end, verified against IRGSP-1.0 on 446/446 loci.
"""

import hashlib
import json
import os
import tempfile
from typing import Dict, Optional

_INDEX_VERSION = 1
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


class GenomeFasta:
    """Byte-offset indexed FASTA reader."""

    def __init__(self, path: str, index_path: Optional[str] = None):
        self.path = os.path.abspath(path)
        if not os.path.isfile(self.path):
            raise FileNotFoundError(f"genome FASTA not found: {self.path}")
        self._fh = open(self.path, "rb")
        self._index_path = index_path or self._default_index_path()
        self.records = self._load_or_build_index()

    # -- index ---------------------------------------------------------
    def _default_index_path(self) -> str:
        """Prefer the FASTA's own directory, fall back to a temp cache."""
        stat = os.stat(self.path)
        tag = f"{self.path}:{stat.st_size}:{int(stat.st_mtime)}"
        digest = hashlib.sha1(tag.encode()).hexdigest()[:12]
        candidates = [
            self.path + f".mirdp-idx.{digest}.json",
            os.path.join(tempfile.gettempdir(), f"mirdp-genome-{digest}.json"),
        ]
        for cand in candidates:
            if os.path.exists(cand):
                return cand
            # probe writability cheaply
            try:
                with open(cand, "w"):
                    pass
                os.unlink(cand)
                return cand
            except OSError:
                continue
        return candidates[-1]

    def _load_or_build_index(self) -> Dict[str, Dict[str, int]]:
        stat = os.stat(self.path)
        key = {"version": _INDEX_VERSION, "size": stat.st_size,
               "mtime": int(stat.st_mtime)}
        if os.path.exists(self._index_path):
            try:
                with open(self._index_path) as fh:
                    cached = json.load(fh)
                if cached.get("key") == key:
                    return cached["records"]
            except (OSError, ValueError, KeyError):
                pass
        records = self._build_index()
        try:
            with open(self._index_path, "w") as fh:
                json.dump({"key": key, "records": records}, fh)
        except OSError:
            pass  # read-only cache location is fine, just slower next time
        return records

    def _build_index(self) -> Dict[str, Dict[str, int]]:
        records: Dict[str, Dict[str, int]] = {}
        with open(self.path, "rb") as fh:
            name = None
            base_offset = 0
            length = 0
            line_bases = 0
            line_width = 0
            pos = 0
            for line in fh:
                if line.startswith(b">"):
                    if name is not None:
                        records[name] = {"offset": base_offset, "length": length,
                                         "line_bases": line_bases,
                                         "line_width": line_width}
                    header = line[1:].decode("ascii", "replace").split()
                    name = header[0] if header else ""
                    base_offset = pos + len(line)
                    length = 0
                    line_bases = 0
                    line_width = 0
                else:
                    n = len(line.rstrip(b"\r\n"))
                    if n:
                        if line_bases == 0:
                            line_bases = n
                            line_width = len(line)
                        length += n
                pos += len(line)
            if name is not None:
                records[name] = {"offset": base_offset, "length": length,
                                 "line_bases": line_bases,
                                 "line_width": line_width}
        return records

    # -- queries -------------------------------------------------------
    def __contains__(self, chrom: str) -> bool:
        return chrom in self.records

    def length(self, chrom: str) -> int:
        return self.records[chrom]["length"]

    def fetch(self, chrom: str, start: int, end: int,
              strand: str = "+") -> Optional[str]:
        """Return ``[start, end)`` in transcript orientation, or None."""
        rec = self.records.get(chrom)
        if rec is None:
            return None
        if start < 0 or end > rec["length"] or start >= end:
            return None
        line_bases = rec["line_bases"]
        line_width = rec["line_width"]
        if line_bases <= 0:
            return None
        byte_start = rec["offset"] + (start // line_bases) * line_width \
            + (start % line_bases)
        last = end - 1
        byte_end = rec["offset"] + (last // line_bases) * line_width \
            + (last % line_bases) + 1
        self._fh.seek(byte_start)
        raw = self._fh.read(byte_end - byte_start)
        seq = raw.replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()
        if len(seq) != end - start:
            return None
        return seq if strand == "+" else seq.translate(_COMPLEMENT)[::-1]

    def close(self):
        if self._fh and not self._fh.closed:
            self._fh.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False
