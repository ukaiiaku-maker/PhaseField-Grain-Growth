#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

from grain_growth_pf.campaign import enumerate_campaign
from grain_growth_pf.config import ModelConfig
from grain_growth_pf.io.provenance import canonical_hash, git_sha
from grain_growth_pf.migration_closure import MigrationClosureSimulation
from grain_growth_pf.pf.initial_conditions import initial_condition_identity, prepare_initial_condition


def _worker(payload: tuple[dict[str, Any], str, str]) -> dict[str, str]:
    config_dict, output, sha = payload
    try:
        simulation = MigrationClosureSimulation(
            ModelConfig.from_dict(config_dict), output, code_sha=sha
        )
        simulation.run()
        return {"path": output, "status": "completed"}
    except BaseException as exc:
        Path(output).mkdir(parents=True, exist_ok=True)
        (Path(output) / "closure_failure.txt").write_text(
            f"{type(exc).__name__}: {exc}\n", encoding="utf-8"
        )
        return {"path": output, "status": "failed", "error": str(exc)}


def load_configs(spec_path: Path, output_root: Path, sha: str) -> list[ModelConfig]:
    spec = yaml.safe_load(spec_path.read_text(encoding="utf-8"))
    configs = enumerate_campaign(spec)
    if not spec.get("prepare_initial_conditions", False):
        return configs

    base = ModelConfig.from_dict(spec.get("base", {}))
    cache_root = output_root.parent / "initial_conditions"
    files: dict[int, str] = {}
    for seed in sorted({config.seed for config in configs}):
        identity = initial_condition_identity(base.pf, seed, base.parameters, sha)
        state = cache_root / f"seed-{seed}-{identity[:16]}.npz"
        prepare_initial_condition(base.pf, seed, base.parameters, state, sha)
        files[seed] = str(state)
    return [replace(
        config,
        parameters={**config.parameters, "initial_state_file": files[config.seed]},
    ) for config in configs]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run matched hybrid/gate-only migration-closure comparisons."
    )
    parser.add_argument("spec")
    parser.add_argument("--output-root", default="results/migration_closure")
    parser.add_argument("--processes", type=int, default=8)
    args = parser.parse_args()

    sha = git_sha()
    output_root = Path(args.output_root)
    configs = load_configs(Path(args.spec), output_root, sha)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    identity = canonical_hash({"configs": [c.to_dict() for c in configs], "sha": sha})[:10]
    root = output_root / f"{stamp}-{identity}"
    root.mkdir(parents=True, exist_ok=False)

    payloads = []
    for config in configs:
        run_dir = root / f"{config.regime}-T{config.pf.temperature:g}-s{config.seed}"
        payloads.append((config.to_dict(), str(run_dir), sha))

    workers = min(max(1, args.processes), len(payloads))
    if workers == 1:
        outcomes = [_worker(payload) for payload in payloads]
    else:
        with mp.get_context("spawn").Pool(workers) as pool:
            outcomes = pool.map(_worker, payloads, chunksize=1)

    status = "completed" if all(item["status"] == "completed" for item in outcomes) else "failed"
    manifest = {
        "status": status,
        "git_sha": sha,
        "source_spec": args.spec,
        "workers": workers,
        "runs": [item["path"] for item in outcomes],
        "outcomes": outcomes,
    }
    (root / "campaign_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(root)
    if status != "completed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
