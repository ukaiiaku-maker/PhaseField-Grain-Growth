#!/usr/bin/env python3
"""Bounded diffuse-PF qualification for the preregistered A2 law."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter

from grain_growth_pf.config import PFConfig
from grain_growth_pf.mechanics.anisotropy import LADDER, angular_normalization
from grain_growth_pf.pf.anisotropic import anisotropic_energy_gradient
from grain_growth_pf.pf.geometry import circular_grain, planar_interface, voronoi_polycrystal
from grain_growth_pf.pf.solver import MultiphaseFieldSolver

CGAMMA = 1.2060485247581734
CMOBILITY = 0.924620319961451


def config(shape, *, dt=0.04, dx=1.0, width=4.0, strength="A2_STRONG"):
    return PFConfig(
        shape=shape, grid_spacing=dx, interface_width=width, time_step=dt,
        gb_energy=1.0, intrinsic_mobility=4.0, adaptive_stepping=True,
        anisotropy_strength=strength,
        anisotropy_energy_normalization=CGAMMA,
        anisotropy_mobility_normalization=CMOBILITY,
        grain_extinction_threshold=0.05,
    )


def evolve_steps(eta, cfg, orientations, steps):
    solver = MultiphaseFieldSolver(eta, cfg, orientations=orientations)
    energies = [solver._anisotropic_energy()] if solver.anisotropic else []
    accepted = []
    for _ in range(steps):
        record = solver.step()
        energies.append(record.interfacial_energy)
        accepted.append(record.dt)
    return solver, np.asarray(energies), np.asarray(accepted)


def refinement_runs(eta, orientations, coarse_steps):
    probe = MultiphaseFieldSolver(
        eta.copy(), config(eta.shape[1:]), orientations=orientations
    )
    stability_limit = probe.stable_dt()
    runs = []
    for factor, steps in ((1.0, coarse_steps), (0.5, 2 * coarse_steps)):
        requested = stability_limit * factor
        solver, energies, accepted = evolve_steps(
            eta.copy(), config(eta.shape[1:], dt=requested), orientations, steps
        )
        runs.append({
            "factor": factor,
            "requested_dt": requested,
            "accepted_dt_min": float(np.min(accepted)),
            "accepted_dt_max": float(np.max(accepted)),
            "accepted_exact": bool(np.all(accepted == requested)),
            "steps": steps,
            "physical_time": solver.time,
            "initial_energy": float(energies[0]),
            "final_energy": float(energies[-1]),
            "maximum_energy_increase": float(np.max(np.diff(energies))),
            "positive_energy_steps": int(np.count_nonzero(np.diff(energies) > 1e-12)),
            "constraint_error": float(np.max(np.abs(solver.eta.sum(axis=0) - 1))),
            "minimum_phase_value": float(np.min(solver.eta)),
            "finite": bool(np.all(np.isfinite(solver.eta))),
            "solver": solver,
        })
    return stability_limit, runs


def serialized_run(run):
    return {key: value for key, value in run.items() if key != "solver"}


def force_gradient_check():
    rng = np.random.default_rng(44)
    eta = rng.uniform(0.1, 1.0, (3, 5, 6))
    eta /= eta.sum(axis=0)
    strength = LADDER["A2_STRONG"]
    args = (
        np.ones(3, dtype=bool), np.array([0.11, 0.57, 1.02]),
        1.0, 4.0, 4.0, 1.0, True, strength.g_min,
        strength.inclination_weight, strength.support_power,
        strength.mobility_exponent, angular_normalization(strength),
        CGAMMA, CMOBILITY, True,
    )
    energy, derivative = anisotropic_energy_gradient(eta, *args)
    epsilon = 1e-7
    errors = []
    for index in ((0, 1, 2), (1, 3, 4), (2, 2, 1)):
        plus = eta.copy(); plus[index] += epsilon
        minus = eta.copy(); minus[index] -= epsilon
        ep = anisotropic_energy_gradient(plus, *args)[0]
        em = anisotropic_energy_gradient(minus, *args)[0]
        fd = (ep - em) / (2 * epsilon)
        errors.append(abs(fd - derivative[index]) / max(abs(fd), 1e-12))
    return {"energy": float(energy), "maximum_relative_error": max(errors)}


def a0_nesting():
    eta = circular_grain((40, 40), 10, 4)
    base_cfg = PFConfig(
        shape=(40, 40), interface_width=4, time_step=0.04,
        intrinsic_mobility=4.0, adaptive_stepping=True,
    )
    a0_cfg = PFConfig(**{
        **base_cfg.__dict__, "anisotropy_strength": "A0_ISOTROPIC",
    })
    base = MultiphaseFieldSolver(eta.copy(), base_cfg)
    nested = MultiphaseFieldSolver(eta.copy(), a0_cfg)
    for _ in range(80):
        left = base.step()
        right = nested.step()
        if left != right or not np.array_equal(base.eta, nested.eta):
            return {"exact": False, "step": base.step_number}
    return {"exact": True, "step": base.step_number, "time": base.time}


def planar_checks():
    rows = []
    for angle in (0.0, np.pi / 8, np.pi / 4, 3 * np.pi / 8):
        eta = planar_interface((48, 48), 4, angle=angle)
        initial_mass = float(eta[1].sum())
        stability_limit, runs = refinement_runs(
            eta, np.array([0.13, 0.71]), coarse_steps=160
        )
        coarse = runs[0]["solver"]
        fine = runs[1]["solver"]
        rows.append({
            "angle": angle,
            "stability_limit": stability_limit,
            "runs": [serialized_run(run) for run in runs],
            "coarse_relative_mass_change": float(
                (coarse.eta[1].sum() - initial_mass) / initial_mass
            ),
            "fine_relative_mass_change": float(
                (fine.eta[1].sum() - initial_mass) / initial_mass
            ),
            "timestep_relative_field_error": float(
                np.linalg.norm(coarse.eta - fine.eta) / np.linalg.norm(fine.eta)
            ),
        })
    return rows


def inclusion_and_refinement():
    eta = circular_grain((48, 48), 12, 4)
    stability_limit, runs = refinement_runs(
        eta, np.array([0.07, 0.83]), coarse_steps=800
    )
    coarse = runs[0]["solver"]
    fine = runs[1]["solver"]
    return {
        "stability_limit": stability_limit,
        "runs": [serialized_run(run) for run in runs],
        "coarse_mass": float(coarse.eta[1].sum()),
        "fine_mass": float(fine.eta[1].sum()),
        "timestep_relative_field_error": float(
            np.linalg.norm(coarse.eta - fine.eta) / np.linalg.norm(fine.eta)
        ),
    }


def triple_junction_check():
    shape = (48, 48)
    yy, xx = np.indices(shape)
    angle = np.mod(np.arctan2(yy - 23.5, xx - 23.5), 2 * np.pi)
    labels = np.floor(3 * angle / (2 * np.pi)).astype(int)
    eta = np.stack([gaussian_filter((labels == i).astype(float), 1.2) for i in range(3)])
    eta /= eta.sum(axis=0, keepdims=True)
    stability_limit, runs = refinement_runs(
        eta, np.array([0.03, 0.49, 1.11]), coarse_steps=1200
    )
    coarse = runs[0]["solver"]
    fine = runs[1]["solver"]
    return {
        "stability_limit": stability_limit,
        "runs": [serialized_run(run) for run in runs],
        "timestep_relative_field_error": float(
            np.linalg.norm(coarse.eta - fine.eta) / np.linalg.norm(fine.eta)
        ),
    }


def restart_and_topology():
    eta, _, orientations = voronoi_polycrystal((40, 40), 12, 93, width=2)
    cfg = config((40, 40))
    continuous = MultiphaseFieldSolver(eta.copy(), cfg, orientations=orientations)
    active_counts = [int(continuous.active_phases.sum())]
    energies = [continuous._anisotropic_energy()]
    records = []
    for _ in range(40):
        record = continuous.step()
        records.append(record)
        energies.append(record.interfacial_energy)
        active_counts.append(int(np.count_nonzero(continuous.active_phases)))
    interrupted = MultiphaseFieldSolver(eta.copy(), cfg, orientations=orientations)
    interrupted.run(17)
    restored = MultiphaseFieldSolver(eta.copy(), cfg, orientations=orientations)
    restored.load_state_dict(interrupted.state_dict())
    restored.run(23)
    return {
        "restart_exact": bool(np.array_equal(continuous.eta, restored.eta)),
        "time_exact": bool(continuous.time == restored.time),
        "initial_active": active_counts[0], "final_active": active_counts[-1],
        "no_resurrection": bool(np.all(np.diff(active_counts) <= 0)),
        "finite": bool(np.all(np.isfinite(continuous.eta))),
        "minimum_phase_value": float(np.min(continuous.eta)),
        "maximum_energy_increase": float(np.max(np.diff(energies))),
        "maximum_constraint_error": max(r.max_constraint_error for r in records),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID") or not os.environ.get("SLURMD_NODENAME"):
        raise RuntimeError("consolidated PF qualification requires an HPC3 allocation")
    args.output.mkdir(parents=True, exist_ok=False)
    report = {
        "schema": "anisotropic-diffuse-pf-consolidated-qualification-v2",
        "job_id": os.environ["SLURM_JOB_ID"],
        "strength": "A2_STRONG", "C_gamma": CGAMMA, "C_M": CMOBILITY,
        "force_gradient": force_gradient_check(),
        "a0_nesting": a0_nesting(),
        "planar": planar_checks(),
        "inclusion": inclusion_and_refinement(),
        "triple_junction": triple_junction_check(),
        "restart_topology": restart_and_topology(),
    }
    planar_runs = [run for row in report["planar"] for run in row["runs"]]
    evolving_runs = (
        planar_runs + report["inclusion"]["runs"]
        + report["triple_junction"]["runs"]
    )
    matched_times = []
    for row in report["planar"]:
        matched_times.append(abs(row["runs"][0]["physical_time"] - row["runs"][1]["physical_time"]) <= 1e-12)
    matched_times.extend((
        abs(report["inclusion"]["runs"][0]["physical_time"] - report["inclusion"]["runs"][1]["physical_time"]) <= 1e-12,
        abs(report["triple_junction"]["runs"][0]["physical_time"] - report["triple_junction"]["runs"][1]["physical_time"]) <= 1e-12,
    ))
    raw_gates = {
        "force_gradient": report["force_gradient"]["maximum_relative_error"] <= 1e-4,
        "a0_exact": report["a0_nesting"]["exact"],
        "two_accepted_timesteps": all(run["accepted_exact"] for run in evolving_runs),
        "matched_physical_times": all(matched_times),
        "manufactured_energy_descent": max(run["maximum_energy_increase"] for run in evolving_runs) <= 1e-9,
        "manufactured_finite": all(run["finite"] for run in evolving_runs),
        "manufactured_nonnegative": min(run["minimum_phase_value"] for run in evolving_runs) >= -1e-12,
        "manufactured_phase_sum": max(run["constraint_error"] for run in evolving_runs) <= 1e-12,
        "restart_exact": report["restart_topology"]["restart_exact"],
        "restart_time_exact": report["restart_topology"]["time_exact"],
        "no_resurrection": report["restart_topology"]["no_resurrection"],
        "topology_energy_descent": report["restart_topology"]["maximum_energy_increase"] <= 1e-9,
        "topology_finite": report["restart_topology"]["finite"],
        "topology_nonnegative": report["restart_topology"]["minimum_phase_value"] >= -1e-12,
        "topology_phase_sum": report["restart_topology"]["maximum_constraint_error"] <= 1e-12,
    }
    gates = {name: bool(value) for name, value in raw_gates.items()}
    report["gates"] = gates
    report["passed"] = all(gates.values())
    (args.output / "qualification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
