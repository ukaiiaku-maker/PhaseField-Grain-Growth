#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from grain_growth_pf.campaign import enumerate_campaign
from grain_growth_pf.config import ModelConfig
from grain_growth_pf.io.provenance import canonical_hash, git_sha
from grain_growth_pf.migration_closure import MigrationClosureSimulation
from grain_growth_pf.pf.initial_conditions import initial_condition_identity, prepare_initial_condition


class ClosureFrameSimulation(MigrationClosureSimulation):
    def _write_frame(self, force: bool = False) -> None:
        cadence = max(1, int(self.config.parameters.get("video_frame_cadence", 10)))
        step = int(self.solver.step_number)
        if not force and step % cadence != 0:
            return
        frame_dir = self.output_dir / "frames"
        frame_dir.mkdir(parents=True, exist_ok=True)
        path = frame_dir / f"frame-{step:07d}.npz"
        if path.exists():
            return
        labels = self.solver.labels.astype(np.uint16, copy=True)
        blocked = np.zeros(self.config.pf.shape, dtype=np.uint8)
        shear = np.zeros(self.config.pf.shape, dtype=np.float32)
        free_volume = np.zeros(self.config.pf.shape, dtype=np.float32)
        shape = np.asarray(self.config.pf.shape)
        for segment in self.snapshot.boundaries.values():
            domain = self.domains.get(segment.entity_id)
            if domain is None or not len(segment.points):
                continue
            points = segment.points.astype(int) % shape
            yy, xx = points[:, 0], points[:, 1]
            if domain.blocked:
                blocked[yy, xx] = 1
            shear[yy, xx] = float(domain.shear.state)
            free_volume[yy, xx] = float(domain.free_volume.deficit)
        np.savez_compressed(
            path, labels=labels, blocked=blocked, shear=shear, free_volume=free_volume,
            mobility=self.solver.mobility_scale.astype(np.float32),
            time=np.asarray(float(self.solver.time)), step=np.asarray(step),
            temperature=np.asarray(float(self.config.pf.temperature)),
        )

    def _save_checkpoint(self) -> None:
        super()._save_checkpoint()
        self._write_frame(force=False)

    def run(self) -> Path:
        self._write_frame(force=True)
        result = super().run()
        self._write_frame(force=True)
        return result


def _worker(payload: tuple[dict[str, Any], str, str]) -> dict[str, str]:
    config_dict, output, sha = payload
    try:
        ClosureFrameSimulation(ModelConfig.from_dict(config_dict), output, code_sha=sha).run()
        return {"path": output, "status": "completed"}
    except BaseException as exc:
        Path(output).mkdir(parents=True, exist_ok=True)
        (Path(output) / "video_failure.txt").write_text(
            f"{type(exc).__name__}: {exc}\n", encoding="utf-8"
        )
        return {"path": output, "status": "failed", "error": str(exc)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "spec", nargs="?", default="configs/production/migration_closure_video_c5.yaml"
    )
    parser.add_argument("--output-root", default="results/migration_closure_video")
    parser.add_argument("--processes", type=int, default=4)
    args = parser.parse_args()

    spec_path = Path(args.spec)
    spec = yaml.safe_load(spec_path.read_text(encoding="utf-8"))
    configs = enumerate_campaign(spec)
    sha = git_sha()
    output_root = Path(args.output_root)

    if spec.get("prepare_initial_conditions", False):
        base = ModelConfig.from_dict(spec.get("base", {}))
        cache_root = output_root.parent / "initial_conditions"
        files: dict[int, str] = {}
        for seed in sorted({config.seed for config in configs}):
            identity = initial_condition_identity(base.pf, seed, base.parameters, sha)
            state = cache_root / f"seed-{seed}-{identity[:16]}.npz"
            prepare_initial_condition(base.pf, seed, base.parameters, state, sha)
            files[seed] = str(state)
        configs = [replace(
            config, parameters={**config.parameters, "initial_state_file": files[config.seed]}
        ) for config in configs]

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    identity = canonical_hash({"configs": [c.to_dict() for c in configs], "sha": sha})[:10]
    root = output_root / f"{stamp}-{identity}"
    root.mkdir(parents=True, exist_ok=False)
    payloads = [
        (config.to_dict(), str(root / f"{config.regime}-T{config.pf.temperature:g}-s{config.seed}"), sha)
        for config in configs
    ]
    workers = min(max(1, args.processes), len(payloads))
    if workers == 1:
        outcomes = [_worker(payload) for payload in payloads]
    else:
        with mp.get_context("spawn").Pool(workers) as pool:
            outcomes = pool.map(_worker, payloads, chunksize=1)
    status = "completed" if all(item["status"] == "completed" for item in outcomes) else "failed"
    (root / "video_manifest.json").write_text(json.dumps({
        "status": status, "git_sha": sha, "source_spec": args.spec,
        "workers": workers, "runs": outcomes,
    }, indent=2) + "\n", encoding="utf-8")
    print(root)
    if status != "completed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
