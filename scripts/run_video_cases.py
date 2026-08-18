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
from grain_growth_pf.pf.initial_conditions import initial_condition_identity, prepare_initial_condition
from grain_growth_pf.simulation import EventResolvedSimulation


class FrameSavingSimulation(EventResolvedSimulation):
    """Production simulation with read-only, restart-independent visualization frames."""

    def _write_video_frame(self, force: bool = False) -> None:
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
            path,
            labels=labels,
            blocked=blocked,
            shear=shear,
            free_volume=free_volume,
            time=np.asarray(float(self.solver.time)),
            step=np.asarray(step),
            temperature=np.asarray(float(self.config.pf.temperature)),
        )

    def _save_checkpoint(self) -> None:
        super()._save_checkpoint()
        self._write_video_frame(force=False)

    def run(self) -> Path:
        self._write_video_frame(force=True)
        result = super().run()
        self._write_video_frame(force=True)
        return result


def _worker(payload: tuple[dict[str, Any], str, str]) -> dict[str, str]:
    config_dict, output, sha = payload
    try:
        sim = FrameSavingSimulation(ModelConfig.from_dict(config_dict), output, code_sha=sha)
        sim.run()
        return {"path": output, "status": "completed"}
    except BaseException as exc:
        Path(output).mkdir(parents=True, exist_ok=True)
        (Path(output) / "video_failure.txt").write_text(
            f"{type(exc).__name__}: {exc}\n", encoding="utf-8"
        )
        return {"path": output, "status": "failed", "error": str(exc)}


def load_configs(spec_path: Path, output_root: Path, sha: str) -> list[ModelConfig]:
    spec = yaml.safe_load(spec_path.read_text(encoding="utf-8"))
    if "regime_catalog" in spec:
        catalog_path = (spec_path.parent / spec["regime_catalog"]).resolve()
        catalog = yaml.safe_load(catalog_path.read_text(encoding="utf-8"))["regimes"]
        names = spec.pop("regime_names", list(catalog))
        spec["regimes"] = {name: catalog[name] for name in names}
    configs = enumerate_campaign(spec)

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
            config,
            parameters={**config.parameters, "initial_state_file": files[config.seed]},
        ) for config in configs]
    return configs


def main() -> None:
    parser = argparse.ArgumentParser(description="Run representative PF cases and preserve movie frames.")
    parser.add_argument("spec", nargs="?", default="configs/production/video_representative_200.yaml")
    parser.add_argument("--output-root", default="results/video_runs")
    parser.add_argument("--processes", type=int, default=6)
    args = parser.parse_args()

    sha = git_sha()
    output_root = Path(args.output_root)
    configs = load_configs(Path(args.spec), output_root, sha)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    identity = canonical_hash({"configs": [config.to_dict() for config in configs], "sha": sha})[:10]
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
        "runs": outcomes,
    }
    (root / "video_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(root)
    if status != "completed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
