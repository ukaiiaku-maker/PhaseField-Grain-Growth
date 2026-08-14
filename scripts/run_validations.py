#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from grain_growth_pf.config import PFConfig
from grain_growth_pf.io.provenance import git_sha, software_versions
from grain_growth_pf.pf.geometry import circular_grain, equivalent_radius, planar_interface
from grain_growth_pf.pf.free_energy import coefficients
from grain_growth_pf.pf.solver import MultiphaseFieldSolver


OUT = Path("results/validation/numerical_validation.json")


def circle_slope(dx: float, mobility: float, physical_time: float = 60.0) -> dict[str, float]:
    size = int(round(64 / dx)); size += size % 2
    cfg = PFConfig(shape=(size, size), grid_spacing=dx, interface_width=4.0,
                   time_step=0.04, intrinsic_mobility=mobility,
                   gb_energy=1.0, adaptive_stepping=True)
    eta = circular_grain(cfg.shape, 14.0 / dx, cfg.interface_width / dx)
    solver = MultiphaseFieldSolver(eta, cfg)
    times, values, energies = [], [], []
    while solver.time < physical_time:
        diag = solver.step()
        if solver.step_number % 30 == 0:
            times.append(solver.time)
            values.append(equivalent_radius(solver.eta[1], dx) ** 2)
            energies.append(diag.interfacial_energy)
    start = max(3, len(times) // 8)
    slope, intercept = np.polyfit(times[start:], values[start:], 1)
    fit = slope * np.asarray(times[start:]) + intercept
    r2 = 1 - np.sum((np.asarray(values[start:]) - fit) ** 2) / np.sum((np.asarray(values[start:]) - np.mean(values[start:])) ** 2)
    return {
        "dx": dx, "mobility": mobility, "slope_R2": float(slope),
        "expected_slope": -2 * mobility, "relative_error": float(abs(slope / (-2 * mobility) - 1)),
        "r_squared": float(r2), "energy_monotonic": bool(np.all(np.diff(energies) <= 1e-9)),
    }


def planar_orientation(angle: float) -> dict[str, float]:
    cfg = PFConfig(shape=(64, 64), grid_spacing=1, interface_width=4,
                   time_step=0.04, intrinsic_mobility=0.2,
                   boundary_conditions="neumann", adaptive_stepping=True)
    initial = planar_interface(cfg.shape, cfg.interface_width, angle)
    solver = MultiphaseFieldSolver(initial.copy(), cfg)
    initial_mass = float(initial[1].sum())
    records = solver.run(400)
    return {
        "angle_degrees": float(np.degrees(angle)),
        "phase_area_change": float(solver.eta[1].sum() - initial_mass),
        "max_field_change": float(np.max(np.abs(solver.eta - initial))),
        "energy_change": float(records[-1].interfacial_energy - records[0].interfacial_energy),
        "constraint_error": records[-1].max_constraint_error,
    }


def planar_surface_energy(angle: float) -> float:
    size, margin, width = 256, 32, 8.0
    eta = planar_interface((size, size), width, angle)
    kappa, well = coefficients(1.0, width)
    gx = np.gradient(eta, axis=2); gy = np.gradient(eta, axis=1)
    density = 0.5 * (0.5 * kappa * (gx * gx + gy * gy) + well * eta**2 * (1 - eta)**2).sum(axis=0)
    side = size - 2 * margin
    contour_length = side / max(abs(np.sin(angle)), abs(np.cos(angle)))
    return float(density[margin:-margin, margin:-margin].sum() / contour_length)


def triple_junction() -> dict[str, object]:
    shape = (72, 72); center = (np.asarray(shape) - 1) / 2
    y, x = np.indices(shape); ry, rx = y - center[0], x - center[1]
    directions = np.array([[np.sin(a), np.cos(a)] for a in (0, 2 * np.pi / 3, 4 * np.pi / 3)])
    scores = np.stack([d[0] * ry + d[1] * rx for d in directions]) / 2.0
    scores -= scores.max(axis=0, keepdims=True)
    eta = np.exp(scores); eta /= eta.sum(axis=0, keepdims=True)
    cfg = PFConfig(shape=shape, interface_width=4, time_step=0.03,
                   intrinsic_mobility=0.15, boundary_conditions="neumann", adaptive_stepping=True)
    solver = MultiphaseFieldSolver(eta, cfg); solver.run(500)
    labels = solver.labels
    rays = []
    for pair in ((0, 1), (0, 2), (1, 2)):
        points = []
        for axis in (0, 1):
            shifted = np.roll(labels, -1, axis=axis)
            mask = ((labels == pair[0]) & (shifted == pair[1])) | ((labels == pair[1]) & (shifted == pair[0]))
            for point in np.argwhere(mask):
                vec = point - center
                radius = np.linalg.norm(vec)
                if 3 < radius < 20:
                    points.append(vec)
        mean = np.asarray(points).mean(axis=0)
        rays.append(float(np.arctan2(mean[0], mean[1]) % (2 * np.pi)))
    rays.sort()
    gaps = np.diff(rays + [rays[0] + 2 * np.pi])
    angles = np.degrees(gaps)
    return {"ray_angles_degrees": np.degrees(rays).tolist(),
            "dihedral_angles_degrees": angles.tolist(),
            "maximum_120_degree_error": float(np.max(np.abs(angles - 120)))}


def main() -> None:
    circles = [circle_slope(dx, mobility) for dx in (1.0, 0.75, 0.5) for mobility in (0.1, 0.2)]
    slopes_at_mobility = {dx: [r["slope_R2"] for r in circles if r["dx"] == dx] for dx in (1.0, 0.75, 0.5)}
    surface_energies = {str(float(np.degrees(a))): planar_surface_energy(a) for a in (0, np.pi / 12, np.pi / 6, np.pi / 4)}
    result = {
        "git_sha": git_sha(), "software": software_versions(),
        "planar_interfaces": [planar_orientation(a) for a in (0, np.pi / 6, np.pi / 4)],
        "planar_surface_energy_by_angle": surface_energies,
        "grid_anisotropy_fraction": float(np.ptp(list(surface_energies.values())) / np.mean(list(surface_energies.values()))),
        "circular_grains": circles,
        "mobility_slope_ratios": {str(dx): values[1] / values[0] for dx, values in slopes_at_mobility.items()},
        "mesh_slope_spread_fraction": float(np.ptp([r["slope_R2"] for r in circles if r["mobility"] == 0.2]) / abs(np.mean([r["slope_R2"] for r in circles if r["mobility"] == 0.2]))),
        "triple_junction": triple_junction(),
    }
    result["passed"] = bool(
        max(r["relative_error"] for r in circles) < 0.12
        and all(r["r_squared"] > 0.999 for r in circles)
        and all(r["energy_monotonic"] for r in circles)
        and max(abs(v - 2) for v in result["mobility_slope_ratios"].values()) < 0.08
        and result["mesh_slope_spread_fraction"] < 0.10
        and result["grid_anisotropy_fraction"] < 0.04
        and result["triple_junction"]["maximum_120_degree_error"] < 8.0
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    if not result["passed"]:
        raise SystemExit("numerical validation failed")


if __name__ == "__main__":
    main()
