#!/usr/bin/env python3
"""Consolidated actual-state and reduced-pilot qualification after correction C."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np

from grain_growth_pf.pf.anisotropic import (
    anisotropic_energy_gradient,
    anisotropic_pairwise_step,
)
from grain_growth_pf.pf.free_energy import free_energy
from grain_growth_pf.pf.geometry import voronoi_polycrystal
from grain_growth_pf.pf.solver import MultiphaseFieldSolver
from scripts.diagnose_anisotropic_first_transition import (
    DT,
    GRAINS,
    SEED,
    SHAPE,
    cfg,
    dense_streaming_oracle,
    energy_args,
    energy_support_entries,
    frozen_energy_gradient,
    global_pair_count,
    law_args,
    support_masks,
)
from scripts.run_anisotropic_reduced_pilot import config, morphology


def field_sha256(values: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(values).view(np.uint8)).hexdigest()


def exact_energy(solver: MultiphaseFieldSolver) -> float:
    if solver.anisotropic:
        return solver._anisotropic_energy()
    c = solver.config
    return float(free_energy(
        solver.eta, c.gb_energy, c.interface_width, c.grid_spacing,
        boundary=c.boundary_conditions,
    ))


def support_summary(eta: np.ndarray, active: np.ndarray) -> dict[str, object]:
    _, support = support_masks(eta, active, True)
    counts = support.sum(axis=0)
    pairs = counts * (counts - 1) // 2
    tails = {}
    for threshold in (1e-4, 1e-6, 1e-8, 1e-10, 1e-12):
        mask = (eta > 0.0) & (eta < threshold)
        tails[f"{threshold:.0e}"] = {
            "point_count": int(np.count_nonzero(mask)),
            "point_mass": float(np.sum(eta[mask])),
        }
    return {
        "active_phase_histogram": np.bincount(counts.ravel(), minlength=GRAINS + 1).tolist(),
        "active_pair_histogram": np.bincount(pairs.ravel()).tolist(),
        "maximum_active_phases": int(np.max(counts)),
        "p95_active_phases": float(np.percentile(counts, 95)),
        "maximum_active_pairs": int(np.max(pairs)),
        "global_active_pair_supports": global_pair_count(support),
        "tails": tails,
    }


def pilot_case(name, eta0, orientations, time_step, steps=64):
    strength = "A0_ISOTROPIC" if name == "A0" else "A2_STRONG"
    solver = MultiphaseFieldSolver(
        eta0.copy(), config(strength, name != "A0", False, time_step),
        orientations=orientations,
    )
    initial_active = solver.active_phases.copy()
    initial_morphology = morphology(eta0)
    energies = [exact_energy(solver)]
    support = [support_summary(solver.eta, solver.active_phases)]
    checkpoint = None
    for index in range(steps):
        solver.step(compute_energy=False)
        energies.append(exact_energy(solver))
        support.append(support_summary(solver.eta, solver.active_phases))
        if index + 1 == steps // 2:
            checkpoint = solver.state_dict()
        if (index + 1) % 16 == 0:
            print(f"{name} {index + 1}/{steps}", flush=True)
    assert checkpoint is not None
    restored = MultiphaseFieldSolver(
        eta0.copy(), config(strength, name != "A0", False, time_step),
        orientations=orientations,
    )
    restored.load_state_dict(checkpoint)
    for _ in range(steps - steps // 2):
        restored.step(compute_energy=False)
    deltas = np.diff(energies)
    tolerance = max(1e-10, abs(energies[0]) * 1e-10)
    final_morphology = morphology(solver.eta)
    validity = {
        "energy_descent": bool(np.max(deltas) <= tolerance),
        "restart_exact": bool(np.array_equal(restored.eta, solver.eta)),
        "restart_time_exact": bool(restored.time == solver.time),
        "finite": bool(np.all(np.isfinite(solver.eta))),
        "nonnegative": bool(np.min(solver.eta) >= -1e-12),
        "phase_sum": bool(np.max(np.abs(solver.eta.sum(axis=0) - 1.0)) <= 1e-12),
        "no_resurrection": bool(not np.any(solver.active_phases & ~initial_active)),
    }
    return {
        "name": name,
        "steps": steps,
        "dt": time_step,
        "physical_time": solver.time,
        "initial_energy": energies[0],
        "final_energy": energies[-1],
        "maximum_energy_increase": float(np.max(deltas)),
        "positive_energy_steps": int(np.count_nonzero(deltas > tolerance)),
        "energy_roundoff_tolerance": tolerance,
        "minimum_phase_value": float(np.min(solver.eta)),
        "maximum_phase_sum_error": float(np.max(np.abs(solver.eta.sum(axis=0) - 1.0))),
        "initial_morphology": initial_morphology,
        "final_morphology": final_morphology,
        "final_field_sha256": field_sha256(solver.eta),
        "restart_field_sha256": field_sha256(restored.eta),
        "support_diagnostics": support,
        "validity": validity,
        "passed": bool(all(validity.values())),
    }


def actual_state_audit(state_path: Path) -> dict[str, object]:
    state = np.load(state_path)
    eta = state["eta_before"]
    active = state["active"]
    orientations = state["orientations"]
    c = cfg()
    energy, gradient = anisotropic_energy_gradient(*energy_args(eta, active, orientations, c))
    oracle_energy, oracle_gradient = dense_streaming_oracle(
        eta, active, orientations, *law_args(c),
    )
    energy_support, _ = support_masks(eta, active, True)
    cell, pi, pj = energy_support_entries(energy_support)
    frozen_energy, frozen_gradient = frozen_energy_gradient(
        eta, orientations, cell, pi, pj, *law_args(c),
    )
    external = np.empty((1, 1, 1), dtype=float)
    result, _, _ = anisotropic_pairwise_step(
        eta, active, orientations, np.ones(SHAPE), external, False, DT,
        *law_args(c), True, False,
    )
    velocity = (result - eta) / DT
    gdot = float(np.sum(gradient * velocity))
    frozen_gdot = float(np.sum(frozen_gradient * velocity))
    directional = []
    for factor in tuple(0.5**power for power in range(11)):
        epsilon = DT * factor
        perturbed = eta + epsilon * velocity
        ep = anisotropic_energy_gradient(*energy_args(perturbed, active, orientations, c))[0]
        ef = frozen_energy_gradient(
            perturbed, orientations, cell, pi, pj, *law_args(c),
        )[0]
        directional.append({
            "epsilon": epsilon,
            "production_quotient": float((ep - energy) / epsilon),
            "frozen_support_quotient": float((ef - frozen_energy) / epsilon),
        })
    trials = []
    for divisor in (1, 2, 4, 8, 16):
        time_step = DT / divisor
        trial, _, _ = anisotropic_pairwise_step(
            eta, active, orientations, np.ones(SHAPE), external, False,
            time_step, *law_args(c), True, False,
        )
        trial_energy = anisotropic_energy_gradient(
            *energy_args(trial, active, orientations, c)
        )[0]
        trials.append({"dt": time_step, "delta_energy": float(trial_energy - energy)})
    gradient_relative_error = float(
        np.linalg.norm(gradient - oracle_gradient)
        / max(np.linalg.norm(oracle_gradient), np.finfo(float).tiny)
    )
    frozen_relative_error = abs(
        directional[-1]["frozen_support_quotient"] - frozen_gdot
    ) / max(1.0, abs(frozen_gdot))
    production_relative_error = abs(
        directional[-1]["production_quotient"] - gdot
    ) / max(1.0, abs(gdot))
    classification = (
        "DERIVATIVE_MATCHES" if production_relative_error <= 1e-3
        else "ACTIVE_BRANCH_NONDIFFERENTIABLE" if frozen_relative_error <= 1e-3
        else "OTHER_IDENTIFIED_DEFECT"
    )
    gates = {
        "dense_oracle_matches": gradient_relative_error <= 1e-10,
        "frozen_branch_derivative_matches": frozen_relative_error <= 1e-3,
        "semidiscrete_direction_descends": gdot < 0.0,
        "all_dt_trials_descend": max(row["delta_energy"] for row in trials) <= 0.0,
    }
    return {
        "source_state": str(state_path),
        "classification": classification,
        "energy": energy,
        "dense_oracle_energy": oracle_energy,
        "gradient_dot_velocity": gdot,
        "frozen_gradient_dot_velocity": frozen_gdot,
        "production_dense_relative_gradient_error": gradient_relative_error,
        "production_directional_relative_error": production_relative_error,
        "frozen_directional_relative_error": frozen_relative_error,
        "directional_derivative": directional,
        "timestep_trials": trials,
        "gates": {key: bool(value) for key, value in gates.items()},
        "passed": bool(all(gates.values())),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("captured_state", type=Path)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("post-fix qualification must run in a Slurm allocation")
    args.output.mkdir(parents=True, exist_ok=False)
    eta0, _, orientations = voronoi_polycrystal(
        SHAPE, GRAINS, SEED, width=2.0, periodic=True,
    )
    report = {
        "schema": "anisotropic-postfix-consolidated-v1",
        "job_id": os.environ["SLURM_JOB_ID"],
        "correction": "POSITIVELY_HOMOGENEOUS_ANISOTROPIC_PAIR_NORM",
        "actual_state": actual_state_audit(args.captured_state),
        "A0": pilot_case("A0", eta0, orientations, DT),
        "A2_energy_only": pilot_case("A2_energy_only", eta0, orientations, DT),
        "A2_energy_only_half_dt": pilot_case(
            "A2_energy_only_half_dt", eta0, orientations, DT / 2.0, steps=128,
        ),
    }
    matched_horizon = abs(
        report["A2_energy_only"]["physical_time"]
        - report["A2_energy_only_half_dt"]["physical_time"]
    ) <= 1e-12
    report["gates"] = {
        "actual_state": report["actual_state"]["passed"],
        "A0": report["A0"]["passed"],
        "A2_energy_only": report["A2_energy_only"]["passed"],
        "A2_half_dt": report["A2_energy_only_half_dt"]["passed"],
        "matched_refinement_horizon": matched_horizon,
    }
    report["passed"] = bool(all(report["gates"].values()))
    (args.output / "postfix-qualification.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps({
        "passed": report["passed"],
        "gates": report["gates"],
        "actual_state_classification": report["actual_state"]["classification"],
    }, indent=2), flush=True)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
