#!/usr/bin/env python3
"""Replay the immutable Phase-1 legacy QIU trajectory with late diagnostics."""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import replace
from pathlib import Path

from grain_growth_pf.config import ModelConfig
from grain_growth_pf.io.checkpoints import atomic_write_text
from grain_growth_pf.io.provenance import git_sha
from run_migration_closure_video import ClosureFrameSimulation


EXPECTED_INITIAL_SHA256 = "106819770289d156b0ff6a35d356a9aa172258f428325ec50c5f62f380ea56f6"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("historical_manifest", type=Path)
    parser.add_argument("initial_state", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--diagnostic-start", type=int, default=9000)
    parser.add_argument("--end-step", type=int, default=10246)
    parser.add_argument("--continue-after-guard", action="store_true")
    args = parser.parse_args()
    if _sha256(args.initial_state) != EXPECTED_INITIAL_SHA256:
        raise SystemExit("initial-state SHA-256 mismatch")
    historical = json.loads(args.historical_manifest.read_text())
    config = ModelConfig.from_dict(historical["config"])
    parameters = {
        **config.parameters,
        "initial_state_file": str(args.initial_state.resolve()),
        "checkpoint_cadence": min(
            500, int(config.parameters.get("checkpoint_cadence", 500))
        ),
    }
    config = replace(
        config, parameters=parameters, max_steps=args.diagnostic_start,
        termination_grains=1,
    )
    # HPC3 invokes this script from a parent staging directory.
    sha = git_sha(Path(__file__).resolve().parents[1])
    args.output.mkdir(parents=True, exist_ok=False)
    run = args.output / "QIU_LEGACY_FORENSIC-T900-s5101"
    ClosureFrameSimulation(config, run, code_sha=sha).run()
    if int(json.loads((run / "checkpoint.json").read_text())["step_number"]) != args.diagnostic_start:
        raise SystemExit("legacy replay did not reach the diagnostic-start checkpoint")

    forensic_parameters = {
        **parameters,
        "qiu_diagnostics_enabled": True,
        "qiu_diagnostic_field_start_step": args.diagnostic_start,
        "qiu_diagnostic_field_cadence": 10,
        "qiu_diagnostic_flush_rows": 16,
        "qiu_guard_terminate": not args.continue_after_guard,
        # The historical output already bounds the onset. These thresholds
        # capture the first unresolved numerical/morphology indicator.
        "qiu_guard_clip_fraction": 0.02,
        "qiu_guard_extinction_count": 10,
        "energy_diagnostic_cadence": 1,
    }
    forensic = replace(
        config, parameters=forensic_parameters, max_steps=args.end_step,
        termination_grains=100,
    )
    ClosureFrameSimulation(forensic, run, resume=True, code_sha=sha).run()
    summary = {
        "schema_version": 1,
        "source_commit": sha,
        "historical_source_commit": historical.get("git_sha"),
        "historical_manifest": str(args.historical_manifest.resolve()),
        "initial_state": str(args.initial_state.resolve()),
        "initial_state_sha256": EXPECTED_INITIAL_SHA256,
        "diagnostic_start_step": args.diagnostic_start,
        "configured_end_step": args.end_step,
        "continue_after_guard": args.continue_after_guard,
        "run_directory": str(run.resolve()),
        "terminal_manifest": json.loads((run / "manifest.json").read_text()),
    }
    atomic_write_text(args.output / "legacy_replay_summary.json", json.dumps(summary, indent=2) + "\n")
    print(args.output.resolve())


if __name__ == "__main__":
    main()
