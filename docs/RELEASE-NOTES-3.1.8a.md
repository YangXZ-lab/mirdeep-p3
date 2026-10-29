# Release notes — miRDeep-P3 3.1.8a

Release date: 2026-11-18
Base: 3.1.7c (commit `8ecfa30`)
Scope: engineering — the identification pipeline becomes a declared list of stages,
every run writes a manifest, and interrupted runs can be resumed.

## 1. What changed

| Ledger | Change |
|---|---|
| O-01 | The 15 inline steps of `identification` are now 12 `Stage` objects in `src/commands/stage_runner.py` |
| O-02 | `identification --resume` skips stages whose recorded artifacts are still valid |
| O-03 | Per-sample `<prefix>.run.json` and batch `<time>.manifest.json` record everything needed to reproduce or audit a run |
| O-08 | `test.sh` gains an analysis smoke test (`Stat` + `Differential_expression`), exit code `6` |
| I-23 | Fixed `analysis DE` / `analysis FA`: the aliases used to print the help and exit 1 |

## 2. Declarative stages

`run_pipeline_for_input()` no longer contains the pipeline body. It builds a
`SampleContext`, then hands a list of stages to `run_stages()`:

```python
stages = build_stages(ctx)          # 12 Stage objects, executed in list order
run_stages(stages, ctx, manifest, resume=args.resume, log_file=log_path)
```

Each `Stage` carries:

| Field | Meaning |
|---|---|
| `sid` | Section number from `docs/MIRDEEP-P3-INTERNALS.md` (e.g. `4.10`) |
| `name` | Short human name (`excise-fold`) |
| `artifacts(ctx)` | Files that must exist for the stage to count as done |
| `run(ctx)` | The work; `raise StageEmpty("reason")` marks an empty (not failed) result |

Adding a step is one list entry; reordering is moving one line; the log line
`[stage] <prefix> 4.10 excise-fold ... ok 49.2s` comes for free.

## 3. `--resume`

```bash
mirdeep-p3 identification -i sample.fq -o out -g genome.fa -t 16 --resume
```

A stage is skipped when all three hold: the manifest records it as `ok`/`skipped`,
the artifacts recorded then equal the artifacts computed now, and those files exist
and are non-empty on disk. Anything else (new input, new parameters, new genome,
new tool set, deleted intermediate) invalidates the manifest and the sample runs from
the start. `--clean` is ignored while resuming and says so.

Cost of the safety check: a resume decision reads one JSON file plus `stat()` calls;
the input md5 is only computed when a manifest is first created.

## 4. Run manifest

`<out>/<prefix>/<prefix>.run.json`:

```json
{
  "schema": 1,
  "version": "3.1.8a",
  "signature": "…",
  "input": {"path": "…", "size": 49000000, "mtime": 1.7e9, "md5": "…", "seq_type": "fastq"},
  "genome": {"path": "…", "size": 31000000, "mtime": 1.7e9},
  "params": {"min_len": 18, "max_len": 26, "rpm_threshold": 5, "max_mappings": 15, "pre_length": 300, "reads_clean": true},
  "tools": {"bowtie": "bowtie-align-s version 1.3.1", "RNAfold": "RNAfold 2.7.2", "samtools": "samtools 1.23.1"},
  "env": {"python": "3.12.13", "platform": "Linux-…"},
  "stages": {"4.10": {"name": "excise-fold", "status": "ok", "artifacts": ["…"], "seconds": 49.24}}
}
```

The batch file `<out>/mirdp3-identification-<time>.manifest.json` lists one entry per
sample (`prefix`, `status`, `message`, `manifest`), so a 100-sample run can be
summarised without opening 100 files.

## 5. Validation

- Mini dataset (`tests/data/mini.fq` + `mini_genome.fasta`): the five key outputs are
  byte-identical to 3.1.7c (`mini_nr_predictions` `3ed4743…`, `.bed` `e17fc65…`,
  `mini_filter_P_prediction` `2e63d35…`, `mini.processed.fa` `760c977…`,
  `mini.total_reads` `2f3610f…`).
- Resume: 12/12 stages skipped on an unchanged run; deleting `mini_nr_predictions.bed`
  re-runs only `4.15` and rebuilds the identical file; editing the signature forces a
  full re-run.
- `./test.sh -t 12 -k`: identification, annotation, target prediction (psrna 120 rows
  deterministic, ssearch 112 rows), `Stat` (4 figures), `Differential_expression`
  (180 miRNAs, 42 DEGs) — `ALL TESTS PASSED`, exit 0.

## 6. Upgrade notes

- No interface removed. `--resume` is opt-in; without it behaviour is unchanged
  (stages still run in the same order, outputs are identical).
- Two new files appear in every run directory. Downstream scripts that glob
  `<out>/*` may need to ignore `*.run.json` / `*.manifest.json`.
- `test.sh` now needs `Rscript` for its last section; without it the section prints
  `[skip]` and the script still exits 0.

## 7. Known limitations

- Batch manifest names have minute resolution, so two runs inside the same minute
  overwrite the batch file (per-sample manifests are unaffected). Tracked as I-24.
- `--resume` trusts size+mtime for the input fingerprint. Touching the file without
  changing it forces a full re-run (safe, just slower); editing it in place while
  keeping size and mtime is not detected (practically impossible by accident).
- Analysis commands are still outside the stage runner (O-04 covered their tool
  resolution, not their structure).

## 8. GitHub sync checklist

1. `git add -A && git commit` in the release tree; tag `v3.1.8a`.
2. Push the commit and the tag to `YangXZ-lab/mirdeep-p3`.
3. If a Docker image is rebuilt, tag it `3.1.8a-full` per `.github/workflows/docker-release.yml`.
4. `CHANGELOG.md` already carries the `[Unreleased] - 3.1.8a` section; move it under a
   dated heading when the tag is cut.
5. Mirror the three docs that changed since 3.1.7c: `docs/MIRDEEP-P3-INTERNALS.md`,
   `docs/MIRDEEP-P3-ISSUES.md`, `docs/RELEASE-NOTES-3.1.8a.md`.
