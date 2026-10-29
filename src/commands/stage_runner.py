#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.
"""
Declarative stage runner for the identification pipeline.

Each stage declares *what it produces* (its artifacts) and *how to produce it*.
Keeping both in one object buys three things the old straight-line function
could not offer:

  * --resume can ask "are this stage's outputs already complete?" instead of
    hard-coding a skip decision per step;
  * the run manifest records per-stage status, timing and artifact sizes, so a
    finished run can be audited without reading the log;
  * adding, removing or reordering a step is a one-line change to the stage list
    in identification.py.

Nothing here knows about miRDeep-P3 specifics: the module is a small, generic
driver (Stage / StageEmpty / RunManifest / run_stages).
"""

import hashlib
import json
import os
import platform
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable, Dict, Iterable, Optional, Sequence

MANIFEST_SCHEMA = 1


class StageEmpty(Exception):
    """A stage legitimately produced no candidate: an empty result, not a failure.

    The reason is written to the log and stored in the manifest; the caller is
    expected to emit the placeholder files downstream steps expect.
    """

    def __init__(self, reason: str):
        super().__init__(reason)
        self.reason = reason


@dataclass
class Stage:
    """One step of the pipeline.

    sid        short identifier matching the numbering in docs (e.g. "4.10")
    name       human-readable name used in logs and the manifest
    artifacts  callable(ctx) -> list of paths the stage must leave behind
    run        callable(ctx) -> None; raises StageEmpty for an empty result
    allow_empty  True when a zero-byte artifact still counts as complete
    """

    sid: str
    name: str
    artifacts: Callable[[object], Sequence[Path]]
    run: Callable[[object], None]
    allow_empty: bool = False


def _now() -> str:
    return datetime.now().astimezone().strftime("%Y-%m-%dT%H:%M:%S%z")


def _append_log(log_file, text: str) -> None:
    if not log_file:
        return
    try:
        with open(log_file, "a") as fh:
            fh.write(text if text.endswith("\n") else text + "\n")
    except OSError:
        pass


def _is_ready(path, allow_empty: bool = False) -> bool:
    try:
        stat = Path(path).stat()
    except OSError:
        return False
    return True if allow_empty else stat.st_size > 0


def artifacts_ready(stage: Stage, ctx) -> bool:
    """True when every declared artifact exists (and is non-empty if required)."""
    artifacts = list(stage.artifacts(ctx))
    if not artifacts:
        return False
    return all(_is_ready(path, stage.allow_empty) for path in artifacts)


def md5_file(path, chunk_size: int = 1 << 20) -> str:
    """MD5 of a file, read in chunks so multi-GB inputs do not blow up memory."""
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk_size), b""):
            digest.update(block)
    return digest.hexdigest()


def digest(payload) -> str:
    """Stable MD5 over a JSON-serialisable payload (used as a run signature)."""
    blob = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    return hashlib.md5(blob).hexdigest()


def file_fingerprint(path) -> Dict:
    """Path plus size/mtime: enough to detect a replaced input without hashing it."""
    target = Path(path)
    try:
        stat = target.stat()
    except OSError:
        return {"path": str(target), "exists": False}
    return {
        "path": str(target),
        "exists": True,
        "size": stat.st_size,
        "mtime": int(stat.st_mtime),
    }


class RunManifest:
    """JSON record of one sample's run: inputs, parameters, tools and stages."""

    def __init__(self, path, data: Dict):
        self.path = Path(path)
        self.data = data
        self.stages: Dict[str, Dict] = data.get("stages") or {}

    # -- construction --------------------------------------------------
    @classmethod
    def create(cls, path, *, tool: str, version: str, step: str, prefix: str,
               input_file, genome, params: Dict, tools: Dict, signature: str,
               input_md5: str = "") -> "RunManifest":
        data = {
            "schema": MANIFEST_SCHEMA,
            "tool": tool,
            "version": version,
            "step": step,
            "prefix": prefix,
            "started_at": _now(),
            "finished_at": None,
            "status": "running",
            "message": "",
            "signature": signature,
            "input": dict(file_fingerprint(input_file),
                          md5=input_md5, seq_type=params.get("type")),
            "genome": file_fingerprint(genome),
            "params": params,
            "tools": tools,
            "env": {
                "python": platform.python_version(),
                "platform": platform.platform(),
            },
            "stages": {},
        }
        manifest = cls(path, data)
        manifest.write()
        return manifest

    @classmethod
    def create_batch(cls, path, *, tool: str, version: str, step: str, params: Dict,
                     tools: Dict, signature: str, inputs: Sequence) -> "RunManifest":
        """Manifest for a whole batch: one entry per sample, no stage list."""
        data = {
            "schema": MANIFEST_SCHEMA,
            "tool": tool,
            "version": version,
            "step": step,
            "prefix": None,
            "started_at": _now(),
            "finished_at": None,
            "status": "running",
            "message": "",
            "signature": signature,
            "inputs": [dict(file_fingerprint(p), md5="", seq_type=params.get("type"))
                       for p in inputs],
            "params": params,
            "tools": tools,
            "env": {
                "python": platform.python_version(),
                "platform": platform.platform(),
            },
            "samples": [],
            "stages": {},
        }
        manifest = cls(path, data)
        manifest.write()
        return manifest

    def add_sample(self, prefix: str, status: str, message: str = "",
                   manifest_path: str = "") -> None:
        self.data.setdefault("samples", []).append({
            "prefix": prefix,
            "status": status,
            "message": message,
            "manifest": manifest_path,
        })
        self.write()

    @classmethod
    def load(cls, path) -> "RunManifest":
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        if not isinstance(data, dict) or "stages" not in data:
            raise ValueError("not a miRDeep-P3 run manifest: %s" % path)
        return cls(path, data)

    # -- queries -------------------------------------------------------
    def signature_matches(self, signature: str) -> bool:
        return self.data.get("signature") == signature

    def is_stage_done(self, stage: Stage, ctx) -> bool:
        """A stage may be skipped only if it finished AND its artifacts survived."""
        record = self.stages.get(stage.sid)
        if not record:
            return False
        if record.get("status") not in ("ok", "skipped"):
            return False
        recorded = [str(Path(p)) for p in (record.get("artifacts") or [])]
        current = [str(Path(p)) for p in stage.artifacts(ctx)]
        if recorded != current:
            # The stage ran with different paths (prefix/parameter change).
            return False
        return artifacts_ready(stage, ctx)

    # -- writes --------------------------------------------------------
    def mark(self, stage: Stage, status: str, artifacts: Iterable,
             message: str = "", seconds: Optional[float] = None) -> None:
        artifacts = [Path(p) for p in artifacts]
        record = {
            "name": stage.name,
            "status": status,
            "artifacts": [str(p) for p in artifacts],
            "artifacts_info": [file_fingerprint(p) for p in artifacts],
            "message": message,
        }
        if seconds is not None:
            record["seconds"] = round(seconds, 2)
        self.stages[stage.sid] = record
        self.data["stages"] = self.stages
        self.write()

    def finish(self, status: str, message: str = "") -> None:
        self.data["status"] = status
        self.data["message"] = message
        self.data["finished_at"] = _now()
        self.write()

    def write(self) -> None:
        """Write atomically: a half-written manifest must never be resumed from."""
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_name(self.path.name + ".part")
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(self.data, fh, indent=2, ensure_ascii=False)
            fh.write("\n")
        os.replace(tmp, self.path)


def run_stages(stages: Sequence[Stage], ctx, manifest: RunManifest,
               resume: bool = False, log_file=None) -> None:
    """Execute the stages in order, recording each outcome in the manifest.

    Raises StageEmpty when a stage reports an empty result; every other exception
    propagates after the failing stage has been recorded as "failed".
    """
    for stage in stages:
        artifacts = [Path(p) for p in stage.artifacts(ctx)]
        if resume and manifest.is_stage_done(stage, ctx):
            _append_log(log_file, "[resume] skip %s %s" % (stage.sid, stage.name))
            manifest.mark(stage, "skipped", artifacts, "already complete")
            print("  [resume] %s: skip %s %s"
                  % (getattr(ctx, "prefix", ""), stage.sid, stage.name), flush=True)
            continue
        _append_log(log_file, "[stage] start %s %s" % (stage.sid, stage.name))
        started = time.time()
        try:
            stage.run(ctx)
        except StageEmpty as exc:
            manifest.mark(stage, "empty", artifacts, exc.reason, time.time() - started)
            _append_log(log_file, "[stage] empty %s %s: %s"
                        % (stage.sid, stage.name, exc.reason))
            raise
        except Exception as exc:                      # noqa: BLE001 - recorded, re-raised
            manifest.mark(stage, "failed", artifacts,
                          "%s: %s" % (type(exc).__name__, exc), time.time() - started)
            _append_log(log_file, "[stage] failed %s %s: %s: %s"
                        % (stage.sid, stage.name, type(exc).__name__, exc))
            raise
        elapsed = time.time() - started
        manifest.mark(stage, "ok", artifacts, "", elapsed)
        _append_log(log_file, "[stage] done %s %s (%.1fs)"
                    % (stage.sid, stage.name, elapsed))
