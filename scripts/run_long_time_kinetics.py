#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import resource
import sys
import time
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from grain_growth_pf.campaign import enumerate_campaign
from grain_growth_pf.config import ModelConfig
from grain_growth_pf.io.checkpoints import atomic_write_text
from grain_growth_pf.io.provenance import canonical_hash, git_sha
from grain_growth_pf.pf.initial_conditions import (
    initial_condition_identity,
    prepare_initial_condition,
)
from run_migration_closure_video import ClosureFrameSimulation


PRIORITY = (
    "B0", "GTSC_GBTJ_Ks025", "GTSC_GBTJ_Ks030_LONG", "GTC_GB",
    "GSC_GBTJ_Ks025", "G", "T", "GT", "QIU",
)


def _directory_bytes(path: Path) -> int:
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def _peak_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


def _worker(payload: tuple[dict[str, Any], str, str, bool]) -> dict[str, Any]:
    raw_config, raw_path, sha, resume = payload
    config = ModelConfig.from_dict(raw_config)
    path = Path(raw_path)
    start = time.perf_counter()
    starting_bytes = _directory_bytes(path) if path.exists() else 0
    try:
        simulation = ClosureFrameSimulation(
            config, path, resume=resume, code_sha=sha
        )
        initialization_seconds = time.perf_counter() - start
        start_step = int(simulation.solver.step_number)
        solver_start = time.perf_counter()
        simulation.run()
        solver_elapsed = time.perf_counter() - solver_start
        elapsed = time.perf_counter() - start
        steps = max(int(simulation.solver.step_number) - start_step, 0)
        size = _directory_bytes(path)
        frame_sizes = [item.stat().st_size for item in (path / "frames").glob("frame-*.npz")]
        trace_size = _directory_bytes(path / "event_traces.parquet") if (path / "event_traces.parquet").exists() else 0
        trace_index_size = _directory_bytes(path / "event_trace_events.parquet") if (path / "event_trace_events.parquet").exists() else 0
        ledger_size = _directory_bytes(path / "events.parquet") if (path / "events.parquet").exists() else 0
        result = {
            "path": str(path), "regime": config.regime, "status": "completed",
            "resumed": resume, "start_step": start_step,
            "end_step": int(simulation.solver.step_number),
            "steps_executed": steps, "wall_seconds": elapsed,
            "initialization_seconds": initialization_seconds,
            "solver_wall_seconds": solver_elapsed,
            "seconds_per_solver_step": solver_elapsed / steps if steps else None,
            "peak_rss_bytes": _peak_rss_bytes(), "run_size_bytes": size,
            "checkpoint_size_bytes": (path / "checkpoint.npz").stat().st_size,
            "ordinary_frame_mean_size_bytes": float(np.mean(frame_sizes)) if frame_sizes else 0.0,
            "ordinary_frame_count": len(frame_sizes),
            "event_trace_size_bytes": trace_size,
            "event_trace_index_size_bytes": trace_index_size,
            "event_ledger_size_bytes": ledger_size,
            "event_trace_growth_bytes_per_step": (trace_size + trace_index_size) / steps if steps else None,
            "bytes_written_this_execution": max(size - starting_bytes, 0),
            "projected_storage_100k_bytes": (
                max(size - starting_bytes, 0) * 100000 / steps if steps else None
            ),
        }
    except BaseException as exc:
        result = {
            "path": str(path), "regime": config.regime, "status": "failed",
            "resumed": resume, "wall_seconds": time.perf_counter() - start,
            "peak_rss_bytes": _peak_rss_bytes(),
            "error": f"{type(exc).__name__}: {exc}",
        }
        path.mkdir(parents=True, exist_ok=True)
        atomic_write_text(path / "long_run_failure.json", json.dumps(result, indent=2) + "\n")
    atomic_write_text(path / "performance.json", json.dumps(result, indent=2) + "\n")
    return result


def _load_new_campaign(args: argparse.Namespace, sha: str) -> tuple[Path, list[ModelConfig], str]:
    spec_path = Path(args.spec)
    spec = yaml.safe_load(spec_path.read_text(encoding="utf-8"))
    configs = enumerate_campaign(spec)
    selected = set(args.regime or (PRIORITY[:2] if args.preflight else PRIORITY))
    unknown = selected.difference(config.regime for config in configs)
    if unknown:
        raise ValueError(f"unknown requested regimes: {sorted(unknown)}")
    configs = [config for name in PRIORITY for config in configs if config.regime == name and name in selected]
    if args.max_steps is not None:
        configs = [replace(config, max_steps=int(args.max_steps)) for config in configs]
    if args.preflight:
        configs = [replace(config, max_steps=int(args.max_steps or 200), termination_grains=1) for config in configs]

    base = ModelConfig.from_dict(spec["base"])
    cache_root = Path(args.output_root).parent / "initial_conditions"
    identity = initial_condition_identity(base.pf, base.seed, base.parameters, sha)
    state = cache_root / f"seed-{base.seed}-{identity[:16]}.npz"
    prepare_initial_condition(base.pf, base.seed, base.parameters, state, sha)
    configs = [replace(
        config, parameters={**config.parameters, "initial_state_file": str(state)}
    ) for config in configs]

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    identity = canonical_hash({"configs": [item.to_dict() for item in configs], "sha": sha})[:10]
    root = Path(args.output_root) / f"{stamp}-{identity}"
    root.mkdir(parents=True, exist_ok=False)
    return root, configs, str(state)


def _load_resume_campaign(root: Path, sha: str) -> tuple[list[ModelConfig], str]:
    manifest = json.loads((root / "campaign_manifest.json").read_text())
    if manifest["git_sha"] != sha:
        raise RuntimeError(
            f"refusing to change source during production: {manifest['git_sha']} != {sha}"
        )
    return [ModelConfig.from_dict(item) for item in manifest["configs"]], manifest["initial_state_file"]


def main() -> None:
    parser = argparse.ArgumentParser(description="Run restart-safe Phase-1 long-time kinetics.")
    parser.add_argument(
        "spec", nargs="?", default="configs/production/long_time_kinetics_900K.yaml"
    )
    parser.add_argument("--output-root", default="results/long_time_kinetics_900K_20260824")
    parser.add_argument("--processes", type=int, default=1)
    parser.add_argument("--regime", action="append")
    parser.add_argument("--max-steps", type=int)
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--resume", type=Path)
    args = parser.parse_args()
    if not 1 <= args.processes <= 2:
        raise SystemExit("Phase 1 permits only one or two concurrent large simulations")

    sha = git_sha()
    if args.resume:
        root = args.resume
        configs, initial_state = _load_resume_campaign(root, sha)
    else:
        root, configs, initial_state = _load_new_campaign(args, sha)

    prior_outcomes: dict[str, dict[str, Any]] = {}
    manifest_path = root / "campaign_manifest.json"
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text())
        prior_outcomes = {item["regime"]: item for item in previous.get("outcomes", [])}

    payloads = []
    retained = []
    for config in configs:
        path = root / f"{config.regime}-T900-s{config.seed}"
        run_manifest = path / "manifest.json"
        completed = False
        if run_manifest.exists():
            data = json.loads(run_manifest.read_text())
            completed = bool(
                data.get("status") == "completed"
                and (
                    int(data.get("steps_completed", -1)) >= config.max_steps
                    or int(data.get("final_grains", config.termination_grains + 1))
                    <= config.termination_grains
                )
            )
        if completed:
            retained.append(prior_outcomes.get(config.regime, {
                "path": str(path), "regime": config.regime,
                "status": "completed", "reused": True,
            }))
            continue
        resume = (path / "checkpoint.npz").exists() and (path / "checkpoint.json").exists()
        payloads.append((config.to_dict(), str(path), sha, resume))

    running = {
        "status": "running", "git_sha": sha, "source_spec": args.spec,
        "initial_state_file": initial_state, "workers": min(args.processes, max(len(payloads), 1)),
        "configs": [config.to_dict() for config in configs],
        "runs": [str(root / f"{config.regime}-T900-s{config.seed}") for config in configs],
        "outcomes": retained,
    }
    atomic_write_text(manifest_path, json.dumps(running, indent=2) + "\n")
    workers = min(args.processes, len(payloads))
    if workers == 0:
        outcomes = []
    elif workers == 1:
        outcomes = [_worker(payload) for payload in payloads]
    else:
        with mp.get_context("spawn").Pool(workers) as pool:
            outcomes = pool.map(_worker, payloads, chunksize=1)
    combined = retained + outcomes
    status = "completed" if all(item["status"] == "completed" for item in combined) else "failed"
    atomic_write_text(manifest_path, json.dumps({
        **running, "status": status, "outcomes": combined,
    }, indent=2) + "\n")
    print(root)
    if status != "completed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
