#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import shutil
import traceback
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

from grain_growth_pf.config import ModelConfig
from grain_growth_pf.io.provenance import canonical_hash, git_sha
from run_migration_closure_video import ClosureFrameSimulation


def _source_runs(root: Path) -> dict[str, Path]:
    manifest_name = "campaign_manifest.json" if (root / "campaign_manifest.json").exists() else "video_manifest.json"
    manifest = json.loads((root / manifest_name).read_text())
    raw_runs = manifest["runs"]
    paths = [Path(item["path"] if isinstance(item, dict) else item) for item in raw_runs]
    result: dict[str, Path] = {}
    for path in paths:
        run = json.loads((path / "manifest.json").read_text())
        regime = str(run["config"]["regime"])
        if regime in result:
            raise ValueError(f"duplicate source regime {regime} under {root}")
        result[regime] = path
    return result


def _worker(payload: tuple[dict[str, Any], str, str, dict[str, Any]]) -> dict[str, Any]:
    config_data, target, sha, restart = payload
    try:
        ClosureFrameSimulation(
            ModelConfig.from_dict(config_data), target, resume=True, code_sha=sha
        ).run()
        manifest_path = Path(target) / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["restart_provenance"] = restart
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        return {"path": target, "status": "completed"}
    except BaseException as exc:
        failure = traceback.format_exc()
        (Path(target) / "extension_failure.txt").write_text(failure)
        return {"path": target, "status": "failed", "error": f"{type(exc).__name__}: {exc}"}


def main() -> None:
    parser = argparse.ArgumentParser(description="Copy and exactly continue corrected closure checkpoints with frames.")
    parser.add_argument("spec")
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--processes", type=int, default=4)
    args = parser.parse_args()

    spec_path = Path(args.spec)
    spec = yaml.safe_load(spec_path.read_text())
    source_root = Path(spec["source_root"])
    sources = _source_runs(source_root)
    mappings = dict(spec["regimes"])
    max_steps = int(spec.get("max_steps", 10000))
    termination_grains = int(spec.get("termination_grains", 70))
    output_cadence = int(spec.get("output_cadence", 10))
    video_cadence = int(spec.get("video_frame_cadence", 50))
    sha = git_sha()
    identity = canonical_hash({
        "source": str(source_root), "regimes": mappings, "max_steps": max_steps,
        "termination_grains": termination_grains, "sha": sha,
    })[:10]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    root = Path(args.output_root) / f"{stamp}-{identity}"
    root.mkdir(parents=True, exist_ok=False)

    payloads = []
    restart_provenance: dict[str, Any] = {}
    for target_regime, source_regime in mappings.items():
        if source_regime not in sources:
            raise ValueError(f"source regime {source_regime!r} is absent from {source_root}")
        source = sources[source_regime]
        source_manifest = json.loads((source / "manifest.json").read_text())
        checkpoint = json.loads((source / "checkpoint.json").read_text())
        config = ModelConfig.from_dict(source_manifest["config"])
        config = replace(
            config,
            regime=str(target_regime),
            max_steps=max_steps,
            termination_grains=termination_grains,
            output_cadence=output_cadence,
            parameters={**config.parameters, "video_frame_cadence": video_cadence},
        )
        target = root / f"{target_regime}-T{config.pf.temperature:g}-s{config.seed}"
        shutil.copytree(source, target)
        restart = {
            "source_run": str(source),
            "source_git_sha": source_manifest["git_sha"],
            "source_config_sha256": source_manifest["config_sha256"],
            "checkpoint_step": int(checkpoint["step_number"]),
            "checkpoint_time": float(checkpoint["time"]),
            "starting_grains": int(source_manifest["final_grains"]),
            "resumed_git_sha": sha,
        }
        restart_provenance[str(target)] = restart
        payloads.append((config.to_dict(), str(target), sha, restart))

    root_manifest = root / "campaign_manifest.json"
    root_manifest.write_text(json.dumps({
        "status": "running", "git_sha": sha, "source_spec": str(spec_path),
        "source_root": str(source_root), "runs": [item[1] for item in payloads],
        "restart_provenance": restart_provenance,
    }, indent=2) + "\n")
    workers = min(max(1, args.processes), len(payloads))
    if workers == 1:
        outcomes = [_worker(payload) for payload in payloads]
    else:
        with mp.get_context("spawn").Pool(workers) as pool:
            outcomes = pool.map(_worker, payloads, chunksize=1)
    status = "completed" if all(item["status"] == "completed" for item in outcomes) else "failed"
    root_manifest.write_text(json.dumps({
        "status": status, "git_sha": sha, "source_spec": str(spec_path),
        "source_root": str(source_root), "workers": workers,
        "runs": [item["path"] for item in outcomes], "outcomes": outcomes,
        "restart_provenance": restart_provenance,
    }, indent=2) + "\n")
    print(root)
    if status != "completed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
