#!/usr/bin/env python3
"""Continue an exact legacy replay checkpoint with the corrected capture guard."""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from grain_growth_pf.config import ModelConfig
from grain_growth_pf.io.checkpoints import atomic_write_text
from grain_growth_pf.io.provenance import git_sha
from run_migration_closure_video import ClosureFrameSimulation


EXPECTED_INITIAL_SHA256 = "106819770289d156b0ff6a35d356a9aa172258f428325ec50c5f62f380ea56f6"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def recovery_step(run: Path) -> int:
    """Validate the paired JSON/NPZ checkpoint and return its solver step."""
    disk_state = json.loads((run / "checkpoint.json").read_text())
    with np.load(run / "checkpoint.npz") as arrays:
        if "checkpoint_state_json" not in arrays:
            raise ValueError("recovery checkpoint lacks embedded atomic metadata")
        embedded_state = json.loads(str(arrays["checkpoint_state_json"]))
    disk_step = int(disk_state["step_number"])
    embedded_step = int(embedded_state["step_number"])
    if disk_step != embedded_step:
        raise ValueError(
            f"checkpoint metadata mismatch: JSON step {disk_step}, NPZ step {embedded_step}"
        )
    return disk_step


def continuation_config(
    historical: dict[str, object], initial_state: Path, *, end_step: int,
) -> ModelConfig:
    """Enable corrected late diagnostics without changing legacy mechanics."""
    config = ModelConfig.from_dict(historical["config"])
    if config.mechanics_backend != "qiu_full_field":
        raise ValueError("legacy continuation requires mechanics_backend=qiu_full_field")
    parameters = {
        **config.parameters,
        "initial_state_file": str(initial_state.resolve()),
        "checkpoint_cadence": min(
            500, int(config.parameters.get("checkpoint_cadence", 500))
        ),
        "qiu_diagnostics_enabled": True,
        "qiu_diagnostic_field_start_step": 9000,
        "qiu_diagnostic_field_cadence": 10,
        "qiu_diagnostic_flush_rows": 16,
        "qiu_guard_terminate": False,
        "qiu_guard_clip_fraction": 0.02,
        "qiu_guard_clip_factor": 1.5,
        "qiu_guard_clip_additive": 0.05,
        "qiu_guard_clip_warmup_steps": 20,
        "qiu_guard_extinction_count": 10,
        "energy_diagnostic_cadence": 1,
    }
    return replace(
        config, parameters=parameters, max_steps=end_step, termination_grains=100,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("historical_manifest", type=Path)
    parser.add_argument("initial_state", type=Path)
    parser.add_argument("run", type=Path)
    parser.add_argument("--expected-step", type=int, default=8000)
    parser.add_argument("--end-step", type=int, default=10246)
    args = parser.parse_args()

    if sha256(args.initial_state) != EXPECTED_INITIAL_SHA256:
        raise SystemExit("initial-state SHA-256 mismatch")
    if not args.run.is_dir():
        raise SystemExit(f"recovered run does not exist: {args.run}")
    actual_step = recovery_step(args.run)
    if actual_step != args.expected_step:
        raise SystemExit(
            f"recovery checkpoint step {actual_step} != expected {args.expected_step}"
        )
    if args.end_step <= actual_step:
        raise SystemExit("end step must be later than the recovery checkpoint")

    historical = json.loads(args.historical_manifest.read_text())
    config = continuation_config(historical, args.initial_state, end_step=args.end_step)
    sha = git_sha(Path(__file__).resolve().parents[1])
    ClosureFrameSimulation(config, args.run, resume=True, code_sha=sha).run()

    terminal = json.loads((args.run / "checkpoint.json").read_text())
    capture_path = args.run / "diagnostic_capture.json"
    capture: dict[str, object] | None = None
    if capture_path.exists():
        capture = json.loads(capture_path.read_text())
        capture_step = int(capture["step"])
        capture.update({
            "scientific_status": "diagnostic_capture_completed_with_post_window",
            "post_capture_terminal_step": int(terminal["step_number"]),
            "post_capture_steps": int(terminal["step_number"]) - capture_step,
        })
        atomic_write_text(capture_path, json.dumps(capture, indent=2) + "\n")

    summary = {
        "schema_version": 1,
        "source_commit": sha,
        "historical_source_commit": historical.get("git_sha"),
        "historical_manifest": str(args.historical_manifest.resolve()),
        "initial_state": str(args.initial_state.resolve()),
        "initial_state_sha256": EXPECTED_INITIAL_SHA256,
        "recovery_step": actual_step,
        "configured_end_step": args.end_step,
        "terminal_step": int(terminal["step_number"]),
        "terminal_time": float(terminal["time"]),
        "diagnostic_capture": capture,
        "run_directory": str(args.run.resolve()),
        "terminal_manifest": json.loads((args.run / "manifest.json").read_text()),
    }
    atomic_write_text(
        args.run.parent / "legacy_replay_continuation_summary.json",
        json.dumps(summary, indent=2) + "\n",
    )
    print(args.run.resolve())


if __name__ == "__main__":
    main()
