#!/usr/bin/env python3
"""Run the preregistered reduced Phase-1 anisotropic polycrystal pilot."""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np

from grain_growth_pf.config import PFConfig
from grain_growth_pf.pf.free_energy import free_energy
from grain_growth_pf.pf.geometry import voronoi_polycrystal
from grain_growth_pf.pf.solver import MultiphaseFieldSolver


SHAPE = (192, 192)
GRAINS = 200
SEED = 5101
STEPS = 128
MIDPOINT = STEPS // 2
ENERGY_NORMALIZATION = 1.2060485247581734
MOBILITY_NORMALIZATION = 0.924620319961451
RESPONSE_THRESHOLD = 0.10
PREREGISTERED_RESPONSES = (
    "phase_field_activity_rms",
    "boundary_density",
    "grain_area_cv",
    "area_weighted_mean_radius",
)


def config(strength: str, energy: bool, mobility: bool, dt: float) -> PFConfig:
    return PFConfig(
        shape=SHAPE,
        grid_spacing=1.0,
        interface_width=4.0,
        time_step=dt,
        gb_energy=1.0,
        intrinsic_mobility=4.0,
        boundary_conditions="periodic",
        adaptive_stepping=False,
        grain_extinction_threshold=0.5,
        anisotropy_strength=strength,
        anisotropic_energy=energy,
        anisotropic_mobility=mobility,
        anisotropy_energy_normalization=ENERGY_NORMALIZATION,
        anisotropy_mobility_normalization=MOBILITY_NORMALIZATION,
    )


def field_sha256(values: np.ndarray) -> str:
    contiguous = np.ascontiguousarray(values)
    return hashlib.sha256(contiguous.view(np.uint8)).hexdigest()


def morphology(eta: np.ndarray) -> dict[str, float | int]:
    labels = np.argmax(eta, axis=0)
    areas = eta.sum(axis=(1, 2))
    present = areas > 1e-10
    present_areas = areas[present]
    boundary_density = 0.5 * float(
        np.mean(labels != np.roll(labels, 1, axis=0))
        + np.mean(labels != np.roll(labels, 1, axis=1))
    )
    isolated = (
        (labels != np.roll(labels, 1, axis=0))
        & (labels != np.roll(labels, -1, axis=0))
        & (labels != np.roll(labels, 1, axis=1))
        & (labels != np.roll(labels, -1, axis=1))
    )
    radii = np.sqrt(present_areas / np.pi)
    return {
        "active_grains": int(np.count_nonzero(present)),
        "boundary_density": boundary_density,
        "isolated_label_fraction": float(np.mean(isolated)),
        "grain_area_cv": float(np.std(present_areas) / np.mean(present_areas)),
        "area_weighted_mean_radius": float(np.sum(present_areas * radii) / np.sum(present_areas)),
    }


def current_energy(solver: MultiphaseFieldSolver) -> float:
    if solver.anisotropic:
        return solver._anisotropic_energy()
    cfg = solver.config
    return float(free_energy(
        solver.eta, cfg.gb_energy, cfg.interface_width, cfg.grid_spacing,
        boundary=cfg.boundary_conditions,
    ))


def advance(
    solver: MultiphaseFieldSolver,
    steps: int,
    energies: list[float] | None = None,
    label: str = "",
) -> None:
    for index in range(steps):
        diagnostic = solver.step(compute_energy=not solver.anisotropic)
        if energies is not None:
            if solver.anisotropic:
                if index:
                    energies.append(float(solver._last_pre_step_energy))
            else:
                energies.append(float(diagnostic.interfacial_energy))
        if (index + 1) % 16 == 0 or index + 1 == steps:
            print(f"{label} {index + 1}/{steps} time={solver.time:.12g}", flush=True)


def run_case(
    name: str,
    eta0: np.ndarray,
    orientations: np.ndarray,
    cfg: PFConfig,
    steps: int,
    output: Path,
) -> dict[str, object]:
    solver = MultiphaseFieldSolver(eta0.copy(), cfg, orientations=orientations)
    initial_active = solver.active_phases.copy()
    initial_morphology = morphology(eta0)
    energies = [current_energy(solver)]

    advance(solver, MIDPOINT if steps == STEPS else steps // 2, energies, name)
    checkpoint = solver.state_dict()
    remaining = steps - solver.step_number
    advance(solver, remaining, energies, name)
    if solver.anisotropic:
        energies.append(current_energy(solver))
    continuous_eta = solver.eta.copy()
    continuous_time = solver.time
    continuous_active = solver.active_phases.copy()

    restored = MultiphaseFieldSolver(eta0.copy(), cfg, orientations=orientations)
    restored.load_state_dict(checkpoint)
    advance(restored, remaining, None, name + "-restart")
    restart_exact = bool(
        np.array_equal(restored.eta, continuous_eta)
        and np.array_equal(restored.active_phases, continuous_active)
        and np.array_equal(restored.orientations, orientations)
    )
    restart_time_exact = bool(restored.time == continuous_time)

    final_morphology = morphology(continuous_eta)
    deltas = np.diff(np.asarray(energies))
    energy_tolerance = max(1e-10, abs(energies[0]) * 1e-10)
    maximum_increase = float(np.max(deltas)) if deltas.size else 0.0
    phase_sum_error = float(np.max(np.abs(continuous_eta.sum(axis=0) - 1.0)))
    minimum = float(np.min(continuous_eta))
    finite = bool(np.all(np.isfinite(continuous_eta)) and np.all(np.isfinite(energies)))
    no_resurrection = bool(not np.any(continuous_active & ~initial_active))
    boundary_proliferation = (
        final_morphology["boundary_density"]
        / max(float(initial_morphology["boundary_density"]), np.finfo(float).tiny)
    )
    isolated_increase = (
        float(final_morphology["isolated_label_fraction"])
        - float(initial_morphology["isolated_label_fraction"])
    )
    stable_morphology = bool(boundary_proliferation <= 1.25 and isolated_increase <= 0.002)
    activity = float(np.sqrt(np.mean(np.sum((continuous_eta - eta0) ** 2, axis=0))))

    with (output / f"{name}-energy.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("state_index", "energy"))
        writer.writerows(enumerate(energies))
    np.savez_compressed(output / f"{name}-labels.npz", labels=np.argmax(continuous_eta, axis=0))

    validity = {
        "finite": finite,
        "nonnegative": minimum >= -1e-12,
        "phase_sum": phase_sum_error <= 1e-12,
        "energy_descent": maximum_increase <= energy_tolerance,
        "restart_exact": restart_exact,
        "restart_time_exact": restart_time_exact,
        "no_resurrection": no_resurrection,
        "stable_morphology": stable_morphology,
    }
    return {
        "name": name,
        "config": asdict(cfg),
        "steps": steps,
        "accepted_dt": cfg.time_step,
        "physical_time": continuous_time,
        "initial_energy": energies[0],
        "final_energy": energies[-1],
        "maximum_energy_increase": maximum_increase,
        "positive_energy_steps": int(np.count_nonzero(deltas > energy_tolerance)),
        "complete_energy_states": len(energies),
        "minimum_phase_value": minimum,
        "maximum_phase_sum_error": phase_sum_error,
        "initial_morphology": initial_morphology,
        "final_morphology": final_morphology,
        "boundary_proliferation_ratio": boundary_proliferation,
        "isolated_label_fraction_increase": isolated_increase,
        "phase_field_activity_rms": activity,
        "final_field_sha256": field_sha256(continuous_eta),
        "restart_field_sha256": field_sha256(restored.eta),
        "validity": validity,
        "passed": bool(all(validity.values())),
    }


def relative_impacts(a0: dict[str, object], combined: dict[str, object]) -> dict[str, float]:
    impacts: dict[str, float] = {}
    for response in PREREGISTERED_RESPONSES:
        if response == "phase_field_activity_rms":
            baseline = float(a0[response])
            value = float(combined[response])
        else:
            baseline = float(a0["final_morphology"][response])
            value = float(combined["final_morphology"][response])
        impacts[response] = abs(value - baseline) / max(abs(baseline), np.finfo(float).tiny)
    return impacts


def main() -> int:
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("the reduced evolving pilot must run in a Slurm allocation")
    output = Path(sys.argv[1] if len(sys.argv) > 1 else "output/reduced-pilot")
    output.mkdir(parents=True, exist_ok=True)
    eta0, seeds, orientations = voronoi_polycrystal(
        SHAPE, GRAINS, SEED, width=2.0, periodic=True
    )
    probe_cfg = config("A2_STRONG", True, True, 1.0)
    probe = MultiphaseFieldSolver(eta0.copy(), probe_cfg, orientations=orientations)
    accepted_dt = probe.stable_dt()

    preregistration = {
        "schema": "anisotropic-reduced-pilot-preregistration-v1",
        "shape": SHAPE,
        "grains": GRAINS,
        "seed": SEED,
        "steps": STEPS,
        "physical_horizon": STEPS * accepted_dt,
        "matched_initial_state": True,
        "matched_accepted_timestep": accepted_dt,
        "cases": ["A0", "A2_energy_only", "A2_mobility_only", "A2_combined"],
        "responses": PREREGISTERED_RESPONSES,
        "consequential_threshold": RESPONSE_THRESHOLD,
        "a3_policy": "run A3 combined exactly once only when all A2 cases are valid and A2 is weak",
        "validity": [
            "finite", "nonnegative", "phase_sum", "energy_descent",
            "restart_exact", "restart_time_exact", "no_resurrection", "stable_morphology",
        ],
    }
    (output / "preregistration.json").write_text(json.dumps(preregistration, indent=2) + "\n")
    np.savez_compressed(output / "matched-initial-state.npz", eta=eta0, seeds=seeds, orientations=orientations)

    specifications = (
        ("A0", "A0_ISOTROPIC", False, False),
        ("A2_energy_only", "A2_STRONG", True, False),
        ("A2_mobility_only", "A2_STRONG", False, True),
        ("A2_combined", "A2_STRONG", True, True),
    )
    cases: dict[str, dict[str, object]] = {}
    for name, strength, energy, mobility in specifications:
        cases[name] = run_case(
            name, eta0, orientations,
            config(strength, energy, mobility, accepted_dt), STEPS, output,
        )
        (output / "progress.json").write_text(json.dumps(cases, indent=2) + "\n")

    impacts = relative_impacts(cases["A0"], cases["A2_combined"])
    all_a2_valid = all(case["passed"] for case in cases.values())
    a2_consequential = max(impacts.values()) >= RESPONSE_THRESHOLD
    a3_evaluated = False
    if all_a2_valid and not a2_consequential:
        a3_evaluated = True
        a3_probe_cfg = config("A3_STRONGER_BOUNDED", True, True, 1.0)
        a3_probe = MultiphaseFieldSolver(eta0.copy(), a3_probe_cfg, orientations=orientations)
        a3_dt = a3_probe.stable_dt()
        horizon = STEPS * accepted_dt
        a3_steps = int(np.ceil(horizon / a3_dt))
        a3_dt = horizon / a3_steps
        cases["A3_combined_once"] = run_case(
            "A3_combined_once", eta0, orientations,
            config("A3_STRONGER_BOUNDED", True, True, a3_dt), a3_steps, output,
        )

    gates = {
        "four_case_matrix_complete": all(name in cases for name, *_ in specifications),
        "all_a2_cases_valid": all_a2_valid,
        "a2_scientifically_consequential": a2_consequential,
        "weak_a2_a3_policy_satisfied": a2_consequential or a3_evaluated,
    }
    report = {
        "schema": "anisotropic-reduced-polycrystal-pilot-v1",
        "job_id": os.environ["SLURM_JOB_ID"],
        "preregistration": preregistration,
        "initial_state_sha256": field_sha256(eta0),
        "cases": cases,
        "combined_a2_relative_impacts": impacts,
        "maximum_combined_a2_relative_impact": max(impacts.values()),
        "a3_evaluated": a3_evaluated,
        "gates": gates,
        "passed": bool(all(gates.values())),
    }
    (output / "pilot.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
