#!/usr/bin/env python3
"""Prepare and hash one additional full-field qualification initial state."""

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
from grain_growth_pf.pf.initial_conditions import (
    initial_condition_identity,
    prepare_initial_condition,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=Path)
    parser.add_argument("seed", type=int)
    parser.add_argument("output_directory", type=Path)
    args = parser.parse_args()
    config = replace(ModelConfig.load(args.config), seed=args.seed)
    if config.regime != "FFT_EIGENSTRAIN_V2":
        raise SystemExit("initial-state preparation requires an FFT_EIGENSTRAIN_V2 config")
    source_sha = git_sha(Path(__file__).resolve().parents[1])
    parameters = dict(config.parameters)
    parameters.pop("initial_state_file", None)
    identity = initial_condition_identity(config.pf, args.seed, parameters, source_sha)
    args.output_directory.mkdir(parents=True, exist_ok=True)
    path = args.output_directory / f"seed-{args.seed}-{identity[:16]}.npz"
    prepare_initial_condition(config.pf, args.seed, parameters, path, source_sha)
    with np.load(path) as data:
        active = int(len(data["orientations"]))
        steps = int(data["equilibration_steps"])
        shape = [active, *config.pf.shape]
    expected = int(parameters.get("equilibrate_to_grains", active))
    if active != expected:
        raise SystemExit(f"prepared {active} active grains, expected exactly {expected}")
    result = {
        "schema_version": 1, "status": "completed",
        "source_commit": source_sha, "config": str(args.config.resolve()),
        "seed": args.seed, "identity": identity,
        "initial_state": str(path.resolve()), "initial_state_sha256": sha256(path),
        "eta_shape": shape, "active_grains": active,
        "equilibration_steps": steps,
        "protocol": {
            key: parameters.get(key) for key in (
                "initial_grains", "equilibration_steps", "equilibrate_to_grains",
                "equilibration_max_steps", "compact_after_equilibration",
            )
        },
    }
    summary = path.with_suffix(".qualification.json")
    atomic_write_text(summary, json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
