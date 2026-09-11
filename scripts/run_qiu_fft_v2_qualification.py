#!/usr/bin/env python3
"""Run or resume one explicitly configured FFT_EIGENSTRAIN_V2 trajectory."""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from grain_growth_pf.config import ModelConfig
from grain_growth_pf.io.checkpoints import atomic_write_text
from grain_growth_pf.io.provenance import git_sha
from grain_growth_pf.simulation import EventResolvedSimulation


class FFTEigenstrainFrameSimulation(EventResolvedSimulation):
    """Base FFT-v2 physics plus compact, restart-safe movie frames."""

    def __init__(self, *args, **kwargs):
        self._last_frame_step = -1
        super().__init__(*args, **kwargs)

    def _write_frame(self, *, force: bool = False) -> None:
        cadence = max(1, int(self.config.parameters.get("video_frame_cadence", 200)))
        if not force and self.solver.step_number % cadence != 0:
            return
        if self._last_frame_step == self.solver.step_number:
            return
        directory = self.output_dir / "frames"; directory.mkdir(exist_ok=True)
        np.savez_compressed(
            directory / f"frame-{self.solver.step_number:07d}.npz",
            labels=self.solver.labels,
            stress=self.full_field.stress,
            eigenstrain=self.full_field.eigenstrain,
            step=np.asarray(self.solver.step_number),
            time=np.asarray(self.solver.time),
            grain_count=np.asarray(len(self.snapshot.grains)),
            G_population=np.asarray(
                np.sqrt(np.prod(self.config.pf.shape) / len(self.snapshot.grains))
            ),
        )
        self._last_frame_step = self.solver.step_number

    def _write_tracks(self) -> None:
        super()._write_tracks()
        self._write_frame()

    def _save_checkpoint(self) -> None:
        super()._save_checkpoint()
        self._write_frame()

    def run(self) -> Path:
        self._write_frame(force=True)
        result = super().run()
        self._write_frame(force=True)
        return result


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=Path)
    parser.add_argument("initial_state", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--max-steps", type=int)
    parser.add_argument("--target", type=float)
    parser.add_argument("--field-cadence", type=int)
    parser.add_argument("--expected-initial-sha256")
    args = parser.parse_args()
    config = ModelConfig.load(args.config)
    if config.regime != "FFT_EIGENSTRAIN_V2":
        raise SystemExit("qualification runner requires the explicit FFT_EIGENSTRAIN_V2 regime")
    initial_hash = sha256(args.initial_state)
    if args.expected_initial_sha256 and initial_hash != args.expected_initial_sha256:
        raise SystemExit("initial-state SHA-256 mismatch")
    parameters = {
        **config.parameters, "initial_state_file": str(args.initial_state.resolve())
    }
    if args.target is not None:
        parameters["external_delta_eta_target"] = args.target
    if args.field_cadence is not None:
        parameters["qiu_diagnostic_field_cadence"] = args.field_cadence
    config = replace(
        config, parameters=parameters,
        max_steps=config.max_steps if args.max_steps is None else args.max_steps,
    )
    resume = (args.output / "checkpoint.npz").exists()
    source_commit = git_sha()
    if resume:
        prior_manifest = json.loads((args.output / "manifest.json").read_text())
        if prior_manifest.get("git_sha") != source_commit:
            raise SystemExit(
                "refusing cross-revision resume: checkpoint was created by "
                f"{prior_manifest.get('git_sha')}, current source is {source_commit}"
            )
    started = time.perf_counter()
    simulation = FFTEigenstrainFrameSimulation(
        config, args.output, resume=resume, code_sha=source_commit
    )
    start_step = simulation.solver.step_number
    simulation.run()
    elapsed = time.perf_counter() - started
    manifest = json.loads((args.output / "manifest.json").read_text())
    usage = resource.getrusage(resource.RUSAGE_SELF)
    peak = int(usage.ru_maxrss)
    summary = {
        "schema_version": 1, "source_commit": source_commit,
        "config": str(args.config.resolve()),
        "initial_state": str(args.initial_state.resolve()),
        "initial_state_sha256": initial_hash,
        "run_directory": str(args.output.resolve()), "resumed": resume,
        "start_step": start_step, "end_step": simulation.solver.step_number,
        "steps_executed": simulation.solver.step_number - start_step,
        "wall_seconds": elapsed,
        "seconds_per_step": elapsed / max(simulation.solver.step_number - start_step, 1),
        "peak_rss_bytes": peak if __import__("sys").platform == "darwin" else peak * 1024,
        "terminal_status": manifest["status"],
        "final_grains": manifest["final_grains"],
        "target": parameters["external_delta_eta_target"],
    }
    atomic_write_text(args.output / "qualification_run_summary.json", json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
