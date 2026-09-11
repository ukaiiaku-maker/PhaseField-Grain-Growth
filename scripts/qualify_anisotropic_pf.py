#!/usr/bin/env python3
"""Bounded diffuse-PF qualification for the preregistered A2 law."""
from __future__ import annotations

import argparse
import json
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


def evolve(eta, cfg, orientations, physical_time):
    solver = MultiphaseFieldSolver(eta, cfg, orientations=orientations)
    energies = [solver._anisotropic_energy()] if solver.anisotropic else []
    steps = 0
    while solver.time < physical_time - 1e-14:
        record = solver.step(dt=min(cfg.time_step, physical_time - solver.time))
        energies.append(record.interfacial_energy)
        steps += 1
    return solver, np.asarray(energies), steps


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
    yy, xx = np.indices((48, 48))
    for angle in (0.0, np.pi / 8, np.pi / 4, 3 * np.pi / 8):
        eta = planar_interface((48, 48), 4, angle=angle)
        initial_mass = float(eta[1].sum())
        solver, energies, steps = evolve(
            eta, config((48, 48)), np.array([0.13, 0.71]), 0.08
        )
        rows.append({
            "angle": angle, "steps": steps,
            "relative_mass_change": float((solver.eta[1].sum() - initial_mass) / initial_mass),
            "maximum_energy_increase": float(np.max(np.diff(energies))),
            "constraint_error": float(np.max(np.abs(solver.eta.sum(axis=0) - 1))),
            "finite": bool(np.all(np.isfinite(solver.eta))),
        })
    return rows


def inclusion_and_refinement():
    eta = circular_grain((48, 48), 12, 4)
    solver, energies, steps = evolve(
        eta, config((48, 48)), np.array([0.07, 0.83]), 0.16
    )
    weights = solver.eta[1]
    yy, xx = np.indices(weights.shape)
    total = weights.sum()
    center = np.array([(yy * weights).sum(), (xx * weights).sum()]) / total
    dy = yy - center[0]; dx = xx - center[1]
    covariance = np.array([
        [(weights * dy * dy).sum(), (weights * dy * dx).sum()],
        [(weights * dy * dx).sum(), (weights * dx * dx).sum()],
    ]) / total
    moments = np.linalg.eigvalsh(covariance)
    aspect = float(np.sqrt(moments[-1] / moments[0]))

    dt_solutions = []
    for requested in (0.04, 0.0004):
        refined, _, refined_steps = evolve(
            eta.copy(), config((48, 48), dt=requested),
            np.array([0.07, 0.83]), 0.04,
        )
        dt_solutions.append((refined.eta, refined_steps))
    dt_error = float(
        np.linalg.norm(dt_solutions[0][0] - dt_solutions[1][0])
        / np.linalg.norm(dt_solutions[1][0])
    )
    return {
        "steps": steps, "initial_energy": float(energies[0]),
        "final_energy": float(energies[-1]),
        "maximum_energy_increase": float(np.max(np.diff(energies))),
        "aspect_ratio": aspect, "timestep_relative_field_error": dt_error,
        "timestep_steps": [item[1] for item in dt_solutions],
    }


def triple_junction_check():
    shape = (48, 48)
    yy, xx = np.indices(shape)
    angle = np.mod(np.arctan2(yy - 23.5, xx - 23.5), 2 * np.pi)
    labels = np.floor(3 * angle / (2 * np.pi)).astype(int)
    eta = np.stack([gaussian_filter((labels == i).astype(float), 1.2) for i in range(3)])
    eta /= eta.sum(axis=0, keepdims=True)
    solver, energies, steps = evolve(
        eta, config(shape), np.array([0.03, 0.49, 1.11]), 0.08
    )
    return {
        "steps": steps, "initial_energy": float(energies[0]),
        "final_energy": float(energies[-1]),
        "maximum_energy_increase": float(np.max(np.diff(energies))),
        "constraint_error": float(np.max(np.abs(solver.eta.sum(axis=0) - 1))),
        "finite": bool(np.all(np.isfinite(solver.eta))),
    }


def restart_and_topology():
    eta, _, orientations = voronoi_polycrystal((40, 40), 12, 93, width=2)
    cfg = config((40, 40))
    continuous = MultiphaseFieldSolver(eta.copy(), cfg, orientations=orientations)
    active_counts = [int(continuous.active_phases.sum())]
    records = continuous.run(40)
    active_counts.extend(int(np.count_nonzero(continuous.active_phases)) for _ in (0,))
    interrupted = MultiphaseFieldSolver(eta.copy(), cfg, orientations=orientations)
    interrupted.run(17)
    restored = MultiphaseFieldSolver(eta.copy(), cfg, orientations=orientations)
    restored.load_state_dict(interrupted.state_dict())
    restored.run(23)
    return {
        "restart_exact": bool(np.array_equal(continuous.eta, restored.eta)),
        "time_exact": bool(continuous.time == restored.time),
        "initial_active": active_counts[0], "final_active": active_counts[-1],
        "no_resurrection": bool(active_counts[-1] <= active_counts[0]),
        "finite": bool(np.all(np.isfinite(continuous.eta))),
        "maximum_constraint_error": max(r.max_constraint_error for r in records),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    report = {
        "schema": "anisotropic-diffuse-pf-qualification-v1",
        "strength": "A2_STRONG", "C_gamma": CGAMMA, "C_M": CMOBILITY,
        "force_gradient": force_gradient_check(),
        "a0_nesting": a0_nesting(),
        "planar": planar_checks(),
        "inclusion": inclusion_and_refinement(),
        "triple_junction": triple_junction_check(),
        "restart_topology": restart_and_topology(),
    }
    gates = {
        "force_gradient": report["force_gradient"]["maximum_relative_error"] <= 1e-4,
        "a0_exact": report["a0_nesting"]["exact"],
        "planar_finite": all(row["finite"] for row in report["planar"]),
        "planar_constraint": max(abs(row["constraint_error"]) for row in report["planar"]) <= 1e-12,
        "planar_energy": max(row["maximum_energy_increase"] for row in report["planar"]) <= 1e-9,
        "inclusion_energy": report["inclusion"]["maximum_energy_increase"] <= 1e-9,
        "timestep_refinement": report["inclusion"]["timestep_relative_field_error"] <= 5e-3,
        "tj_energy": report["triple_junction"]["maximum_energy_increase"] <= 1e-9,
        "tj_finite": report["triple_junction"]["finite"],
        "restart_exact": report["restart_topology"]["restart_exact"],
        "no_resurrection": report["restart_topology"]["no_resurrection"],
    }
    report["gates"] = gates
    report["passed"] = all(gates.values())
    (args.output / "qualification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
