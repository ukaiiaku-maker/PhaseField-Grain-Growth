#!/usr/bin/env python3
"""HPC3-only timestep refinement of the failed diffuse A2 TJ energy gate."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter

from grain_growth_pf.config import PFConfig
from grain_growth_pf.pf.solver import MultiphaseFieldSolver


def initial_tj(shape=(48, 48)):
    yy, xx = np.indices(shape)
    angle = np.mod(np.arctan2(yy - 23.5, xx - 23.5), 2 * np.pi)
    labels = np.floor(3 * angle / (2 * np.pi)).astype(int)
    eta = np.stack([
        gaussian_filter((labels == phase).astype(float), 1.2)
        for phase in range(3)
    ])
    eta /= eta.sum(axis=0, keepdims=True)
    return eta


def make_config(dt):
    return PFConfig(
        shape=(48, 48), interface_width=4, time_step=dt,
        gb_energy=1.0, intrinsic_mobility=4.0, adaptive_stepping=True,
        anisotropy_strength="A2_STRONG",
        anisotropy_energy_normalization=1.2060485247581734,
        anisotropy_mobility_normalization=0.924620319961451,
        grain_extinction_threshold=0.05,
    )


def run(factor, horizon):
    eta = initial_tj()
    orientations = np.array([0.03, 0.49, 1.11])
    probe = MultiphaseFieldSolver(eta.copy(), make_config(0.04), orientations=orientations)
    stability_limit = probe.stable_dt()
    requested = stability_limit * factor
    solver = MultiphaseFieldSolver(
        eta.copy(), make_config(requested), orientations=orientations
    )
    previous = solver._anisotropic_energy()
    initial_energy = previous
    maximum_increase = -np.inf
    positive_steps = 0
    largest_step = None
    while solver.time < horizon - 1e-14:
        record = solver.step(dt=min(requested, horizon - solver.time))
        increase = record.interfacial_energy - previous
        if increase > maximum_increase:
            maximum_increase = increase
            largest_step = {
                "step": record.step, "time": record.time,
                "dt": record.dt, "energy_before": previous,
                "energy_after": record.interfacial_energy,
            }
        if increase > 1e-12:
            positive_steps += 1
        previous = record.interfacial_energy
    return {
        "factor": factor, "stability_limit": stability_limit,
        "requested_dt": requested, "steps": solver.step_number,
        "initial_energy": initial_energy, "final_energy": previous,
        "maximum_energy_increase": maximum_increase,
        "positive_energy_steps": positive_steps,
        "largest_increase_step": largest_step,
        "constraint_error": float(np.max(np.abs(solver.eta.sum(axis=0) - 1))),
        "finite": bool(np.all(np.isfinite(solver.eta))),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--horizon", type=float, default=0.08)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID") or not os.environ.get("SLURMD_NODENAME"):
        raise RuntimeError("TJ refinement requires an HPC3 compute allocation")
    args.output.mkdir(parents=True, exist_ok=False)
    rows = [run(factor, args.horizon) for factor in (1.0, 0.5, 0.25, 0.125)]
    positive = np.asarray([max(row["maximum_energy_increase"], 0.0) for row in rows])
    report = {
        "schema": "anisotropic-diffuse-tj-timestep-refinement-v1",
        "job_id": os.environ["SLURM_JOB_ID"], "rows": rows,
        "monotone_at_finest": bool(rows[-1]["maximum_energy_increase"] <= 1e-9),
        "positive_spike_decreases_with_dt": bool(np.all(np.diff(positive) <= 1e-12)),
    }
    report["passed"] = bool(
        report["monotone_at_finest"] and report["positive_spike_decreases_with_dt"]
    )
    (args.output / "tj-dt-refinement.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
