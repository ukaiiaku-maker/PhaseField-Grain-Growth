#!/usr/bin/env python3
"""Capture and audit the first positive A2 energy-only pilot transition."""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from numba import njit

from grain_growth_pf.config import PFConfig
from grain_growth_pf.mechanics.anisotropy import LADDER, angular_normalization
from grain_growth_pf.pf.anisotropic import (
    _pair_law,
    anisotropic_energy_components,
    anisotropic_energy_gradient,
    anisotropic_pairwise_audit,
    anisotropic_pairwise_step,
)
from grain_growth_pf.pf.geometry import voronoi_polycrystal
from grain_growth_pf.pf.solver import MultiphaseFieldSolver


SHAPE = (192, 192)
GRAINS = 200
SEED = 5101
DT = 6.36493341483619e-5
ENERGY_NORMALIZATION = 1.2060485247581734
MOBILITY_NORMALIZATION = 0.924620319961451
TAIL_THRESHOLDS = (1e-4, 1e-6, 1e-8, 1e-10, 1e-12)
EPSILON_FACTORS = tuple(0.5**power for power in range(11))


def cfg() -> PFConfig:
    return PFConfig(
        shape=SHAPE, grid_spacing=1.0, interface_width=4.0,
        time_step=DT, gb_energy=1.0, intrinsic_mobility=4.0,
        boundary_conditions="periodic", adaptive_stepping=False,
        grain_extinction_threshold=0.5, anisotropy_strength="A2_STRONG",
        anisotropic_energy=True, anisotropic_mobility=False,
        anisotropy_energy_normalization=ENERGY_NORMALIZATION,
        anisotropy_mobility_normalization=MOBILITY_NORMALIZATION,
    )


def law_args(config: PFConfig) -> tuple[object, ...]:
    strength = LADDER["A2_STRONG"]
    return (
        config.gb_energy, config.intrinsic_mobility, config.interface_width,
        config.grid_spacing, True, strength.g_min, strength.inclination_weight,
        strength.support_power, strength.mobility_exponent,
        angular_normalization(strength), config.anisotropy_energy_normalization,
        config.anisotropy_mobility_normalization,
    )


def energy_args(eta, active, orientations, config):
    return (eta, active, orientations, *law_args(config), True)


@njit(cache=True)
def support_masks(eta, active, periodic):
    phases, height, width = eta.shape
    update = np.zeros_like(eta, dtype=np.bool_)
    energy = np.zeros_like(eta, dtype=np.bool_)
    tolerance = 1e-14
    for phase in range(phases):
        if not active[phase]:
            continue
        for y in range(height):
            ym = (y - 1) % height if periodic else max(y - 1, 0)
            yp = (y + 1) % height if periodic else min(y + 1, height - 1)
            for x in range(width):
                xm = (x - 1) % width if periodic else max(x - 1, 0)
                xp = (x + 1) % width if periodic else min(x + 1, width - 1)
                energy[phase, y, x] = (
                    eta[phase, y, x] > tolerance
                    or eta[phase, yp, x] > tolerance
                    or eta[phase, y, xp] > tolerance
                )
                update[phase, y, x] = (
                    energy[phase, y, x]
                    or eta[phase, ym, x] > tolerance
                    or eta[phase, y, xm] > tolerance
                )
    return energy, update


@njit(cache=True)
def global_pair_matrix(support):
    phases, height, width = support.shape
    pairs = np.zeros((phases, phases), dtype=np.bool_)
    local = np.empty(phases, dtype=np.int64)
    for y in range(height):
        for x in range(width):
            count = 0
            for phase in range(phases):
                if support[phase, y, x]:
                    local[count] = phase
                    count += 1
            for left in range(count - 1):
                i = local[left]
                for right in range(left + 1, count):
                    pairs[i, local[right]] = True
    return pairs


def global_pair_count(support):
    return int(np.count_nonzero(global_pair_matrix(support)))


@njit(cache=True)
def energy_support_entries(support):
    phases, height, width = support.shape
    total = 0
    counts = np.sum(support, axis=0)
    for y in range(height):
        for x in range(width):
            total += counts[y, x] * (counts[y, x] - 1) // 2
    cell = np.empty(total, dtype=np.int32)
    phase_i = np.empty(total, dtype=np.int16)
    phase_j = np.empty(total, dtype=np.int16)
    local = np.empty(phases, dtype=np.int64)
    entry = 0
    for y in range(height):
        for x in range(width):
            count = 0
            for phase in range(phases):
                if support[phase, y, x]:
                    local[count] = phase
                    count += 1
            for left in range(count - 1):
                for right in range(left + 1, count):
                    cell[entry] = y * width + x
                    phase_i[entry] = local[left]
                    phase_j[entry] = local[right]
                    entry += 1
    return cell, phase_i, phase_j


@njit(cache=True)
def frozen_energy_gradient(
    eta, orientations, cell, phase_i, phase_j, gamma0, mobility0, width_value,
    dx, periodic, g_min, inclination_weight, support_power,
    mobility_exponent, angular_scale, energy_normalization,
    mobility_normalization,
):
    phases, height, width = eta.shape
    derivative = np.zeros_like(eta)
    energy = 0.0
    gradient_scale = width_value * width_value / (np.pi * np.pi)
    density_scale = 4.0 / width_value
    tolerance = 1e-14
    for entry in range(len(cell)):
        y = cell[entry] // width
        x = cell[entry] - y * width
        yp = (y + 1) % height if periodic else min(y + 1, height - 1)
        xp = (x + 1) % width if periodic else min(x + 1, width - 1)
        i = int(phase_i[entry])
        j = int(phase_j[entry])
        ui = eta[i, y, x]
        uj = eta[j, y, x]
        gix = (eta[i, y, xp] - ui) / dx
        giy = (eta[i, yp, x] - ui) / dx
        gjx = (eta[j, y, xp] - uj) / dx
        gjy = (eta[j, yp, x] - uj) / dx
        px = gix - gjx
        py = giy - gjy
        magnitude = np.sqrt(px * px + py * py)
        theta = np.arctan2(py, px) if magnitude > tolerance else 0.0
        gamma, gamma_first, _ = _pair_law(
            theta, orientations[i], orientations[j], gamma0, mobility0,
            g_min, inclination_weight, support_power, mobility_exponent,
            angular_scale, energy_normalization, mobility_normalization,
        )
        cross = gix * gjx + giy * gjy
        base = ui * uj - gradient_scale * cross
        energy += density_scale * gamma * base * dx * dx
        dgamma_x = 0.0
        dgamma_y = 0.0
        if magnitude > tolerance:
            dgamma_x = -gamma_first * py / (magnitude * magnitude)
            dgamma_y = gamma_first * px / (magnitude * magnitude)
        flux_ix = density_scale * (-gamma * gradient_scale * gjx + base * dgamma_x)
        flux_iy = density_scale * (-gamma * gradient_scale * gjy + base * dgamma_y)
        flux_jx = density_scale * (-gamma * gradient_scale * gix - base * dgamma_x)
        flux_jy = density_scale * (-gamma * gradient_scale * giy - base * dgamma_y)
        derivative[i, y, x] += density_scale * gamma * uj * dx * dx
        derivative[j, y, x] += density_scale * gamma * ui * dx * dx
        derivative[i, y, x] -= (flux_ix + flux_iy) * dx
        derivative[i, y, xp] += flux_ix * dx
        derivative[i, yp, x] += flux_iy * dx
        derivative[j, y, x] -= (flux_jx + flux_jy) * dx
        derivative[j, y, xp] += flux_jx * dx
        derivative[j, yp, x] += flux_jy * dx
    return energy, derivative


@njit(cache=True)
def dense_streaming_oracle(
    eta, active, orientations, gamma0, mobility0, width_value, dx, periodic,
    g_min, inclination_weight, support_power, mobility_exponent, angular_scale,
    energy_normalization, mobility_normalization,
):
    """Independent all-global-pair oracle with no production support list."""
    phases, height, width = eta.shape
    derivative = np.zeros_like(eta)
    total = 0.0
    gradient_scale = width_value * width_value / (np.pi * np.pi)
    density_scale = 4.0 / width_value
    tolerance = 1e-14
    for i in range(phases - 1):
        if not active[i]:
            continue
        for j in range(i + 1, phases):
            if not active[j]:
                continue
            for y in range(height):
                yp = (y + 1) % height if periodic else min(y + 1, height - 1)
                for x in range(width):
                    xp = (x + 1) % width if periodic else min(x + 1, width - 1)
                    if not (
                        (eta[i, y, x] > tolerance or eta[i, yp, x] > tolerance or eta[i, y, xp] > tolerance)
                        and (eta[j, y, x] > tolerance or eta[j, yp, x] > tolerance or eta[j, y, xp] > tolerance)
                    ):
                        continue
                    ui = eta[i, y, x]
                    uj = eta[j, y, x]
                    gix = (eta[i, y, xp] - ui) / dx
                    giy = (eta[i, yp, x] - ui) / dx
                    gjx = (eta[j, y, xp] - uj) / dx
                    gjy = (eta[j, yp, x] - uj) / dx
                    px = gix - gjx
                    py = giy - gjy
                    magnitude = np.sqrt(px * px + py * py)
                    theta = np.arctan2(py, px) if magnitude > tolerance else 0.0
                    gamma, gamma_first, _ = _pair_law(
                        theta, orientations[i], orientations[j], gamma0, mobility0,
                        g_min, inclination_weight, support_power, mobility_exponent,
                        angular_scale, energy_normalization, mobility_normalization,
                    )
                    cross = gix * gjx + giy * gjy
                    base = ui * uj - gradient_scale * cross
                    total += density_scale * gamma * base * dx * dx
                    dgamma_x = 0.0
                    dgamma_y = 0.0
                    if magnitude > tolerance:
                        dgamma_x = -gamma_first * py / (magnitude * magnitude)
                        dgamma_y = gamma_first * px / (magnitude * magnitude)
                    flux_ix = density_scale * (-gamma * gradient_scale * gjx + base * dgamma_x)
                    flux_iy = density_scale * (-gamma * gradient_scale * gjy + base * dgamma_y)
                    flux_jx = density_scale * (-gamma * gradient_scale * gix - base * dgamma_x)
                    flux_jy = density_scale * (-gamma * gradient_scale * giy - base * dgamma_y)
                    derivative[i, y, x] += density_scale * gamma * uj * dx * dx
                    derivative[j, y, x] += density_scale * gamma * ui * dx * dx
                    derivative[i, y, x] -= (flux_ix + flux_iy) * dx
                    derivative[i, y, xp] += flux_ix * dx
                    derivative[i, yp, x] += flux_iy * dx
                    derivative[j, y, x] -= (flux_jx + flux_jy) * dx
                    derivative[j, y, xp] += flux_jx * dx
                    derivative[j, yp, x] += flux_jy * dx
    return total, derivative


def support_statistics(eta, active):
    energy_support, update_support = support_masks(eta, active, True)
    counts = update_support.sum(axis=0)
    pairs = counts * (counts - 1) // 2
    phase_mass = eta.sum(axis=(1, 2))
    tails = {}
    for threshold in TAIL_THRESHOLDS:
        small_values = (eta > 0.0) & (eta < threshold)
        tails[f"{threshold:.0e}"] = {
            "point_count": int(np.count_nonzero(small_values)),
            "point_mass": float(np.sum(eta[small_values])),
            "phase_count_below_total_mass": int(np.count_nonzero((phase_mass > 0.0) & (phase_mass < threshold))),
            "phase_total_mass_below": float(np.sum(phase_mass[(phase_mass > 0.0) & (phase_mass < threshold)])),
        }
    return {
        "active_phase_histogram": np.bincount(counts.ravel(), minlength=GRAINS + 1).tolist(),
        "active_pair_histogram": np.bincount(pairs.ravel()).tolist(),
        "global_active_pair_supports": int(global_pair_count(update_support)),
        "maximum_active_phases": int(np.max(counts)),
        "p95_active_phases": float(np.percentile(counts, 95)),
        "maximum_active_pairs": int(np.max(pairs)),
        "tails": tails,
    }, energy_support, update_support


def frozen_update(eta, chemical, cell, phase_i, phase_j, dt, config):
    height, width = eta.shape[1:]
    flat_cell = cell.astype(np.int64)
    y = flat_cell // width
    x = flat_cell % width
    i = phase_i.astype(np.int64)
    j = phase_j.astype(np.int64)
    difference = chemical[j, y, x] - chemical[i, y, x]
    kinetic = np.pi**2 / (4.0 * config.interface_width * config.grid_spacing**2)
    raw = config.intrinsic_mobility * kinetic * difference
    outgoing = np.zeros_like(eta)
    positive = raw > 0.0
    np.add.at(outgoing, (j[positive], y[positive], x[positive]), raw[positive])
    np.add.at(outgoing, (i[~positive], y[~positive], x[~positive]), -raw[~positive])
    donor = np.ones_like(eta)
    demand = dt * outgoing
    limited_nodes = demand > eta
    donor[limited_nodes] = eta[limited_nodes] / demand[limited_nodes]
    scale = np.where(positive, donor[j, y, x], donor[i, y, x])
    exchange = raw * scale
    rate = np.zeros_like(eta)
    np.add.at(rate, (i, y, x), exchange)
    np.add.at(rate, (j, y, x), -exchange)
    result = eta + dt * rate
    result[np.abs(result) <= 1e-14] = 0.0
    residual = 1.0 - result.sum(axis=0)
    largest = np.argmax(result, axis=0)
    yy, xx = np.indices((height, width))
    result[largest, yy, xx] += residual
    return result, rate, raw, exchange, donor


def graph_metrics(cell, phase_i, phase_j, difference, raw, limited, coefficient, local_count, donor):
    phases = donor.shape[0]
    pixels = donor.shape[1] * donor.shape[2]
    degree = np.zeros(phases * pixels, dtype=float)
    i_index = phase_i.astype(np.int64) * pixels + cell.astype(np.int64)
    j_index = phase_j.astype(np.int64) * pixels + cell.astype(np.int64)
    np.add.at(degree, i_index, coefficient)
    np.add.at(degree, j_index, coefficient)
    limited_mask = limited != raw
    nonzero_gradient = np.abs(difference[np.nonzero(difference)])
    return {
        "maximum_local_active_phase_count": int(np.max(local_count)),
        "p95_local_active_phase_count": float(np.percentile(local_count, 95)),
        "maximum_local_pair_count": int(np.max(local_count * (local_count - 1) // 2)),
        "maximum_weighted_graph_degree": float(np.max(degree)),
        "gershgorin_bound": float(2.0 * np.max(degree)),
        "maximum_mu_difference": float(np.max(np.abs(difference))),
        "p99_mu_difference": float(np.percentile(np.abs(difference), 99)),
        "fraction_donor_limited_exchanges": float(np.mean(limited_mask)),
        "maximum_donor_scaling": float(np.max(donor)),
        "minimum_donor_factor": float(np.min(donor)),
        "maximum_donor_reduction": float(1.0 - np.min(donor)),
        "minimum_nonzero_mu_difference": float(np.min(nonzero_gradient)) if nonzero_gradient.size else 0.0,
        "coefficient_minimum": float(np.min(coefficient)),
        "pair_representation_maximum_error": float(np.max(np.abs(limited - coefficient * difference))),
    }


def trial(pre, active, orientations, dt, config, base_supports):
    external = np.empty((1, 1, 1), dtype=float)
    result, _, _ = anisotropic_pairwise_step(
        pre, active, orientations, np.ones(SHAPE), external, False, dt,
        *law_args(config), True, False,
    )
    audit = anisotropic_pairwise_audit(
        pre, active, orientations, np.ones(SHAPE), dt,
        *law_args(config), True, False,
    )
    chemical, rate, donor, local_count, cell, pi, pj, difference, raw, limited, coefficient, magnitude = audit
    before_energy = anisotropic_energy_gradient(*energy_args(pre, active, orientations, config))[0]
    after_energy = anisotropic_energy_gradient(*energy_args(result, active, orientations, config))[0]
    executed_velocity = (result - pre) / dt
    predicted = -float(np.sum(coefficient * difference * difference))
    metrics = graph_metrics(cell, pi, pj, difference, raw, limited, coefficient, local_count, donor)
    metrics.update({
        "dt": dt,
        "delta_energy": float(after_energy - before_energy),
        "delta_energy_over_dt": float((after_energy - before_energy) / dt),
        "gradient_dot_executed_velocity": float(np.sum(chemical * executed_velocity)),
        "gradient_dot_pair_rate": float(np.sum(chemical * rate)),
        "predicted_pair_dissipation": predicted,
        "minimum_nonzero_pair_gradient": float(np.min(magnitude[magnitude > 0.0])) if np.any(magnitude > 0.0) else 0.0,
        "inclination_below_1e-14": int(np.count_nonzero(magnitude <= 1e-14)),
        "inclination_below_1e-12": int(np.count_nonzero(magnitude <= 1e-12)),
        "inclination_below_1e-10": int(np.count_nonzero(magnitude <= 1e-10)),
        "inclination_below_update_angle_tolerance": int(np.count_nonzero(magnitude * magnitude <= 1e-14)),
    })
    _, after_update_support = support_masks(result, active, True)
    before_energy_support, before_update_support = base_supports
    after_energy_support, _ = support_masks(result, active, True)
    metrics["phase_support_additions"] = int(np.count_nonzero(after_update_support & ~before_update_support))
    metrics["phase_support_removals"] = int(np.count_nonzero(before_update_support & ~after_update_support))
    before_pairs = global_pair_matrix(before_update_support)
    after_pairs = global_pair_matrix(after_update_support)
    metrics["global_pair_support_additions"] = int(np.count_nonzero(after_pairs & ~before_pairs))
    metrics["global_pair_support_removals"] = int(np.count_nonzero(before_pairs & ~after_pairs))

    ecell, ei, ej = energy_support_entries(before_energy_support)
    frozen_before, frozen_gradient = frozen_energy_gradient(
        pre, orientations, ecell, ei, ej, *law_args(config),
    )
    frozen_result, frozen_rate, _, frozen_limited, frozen_donor = frozen_update(
        pre, frozen_gradient, cell, pi, pj, dt, config
    )
    frozen_after = frozen_energy_gradient(
        frozen_result, orientations, ecell, ei, ej, *law_args(config),
    )[0]
    frozen_difference = frozen_gradient[pj.astype(np.int64), cell // SHAPE[1], cell % SHAPE[1]] - frozen_gradient[pi.astype(np.int64), cell // SHAPE[1], cell % SHAPE[1]]
    frozen_coefficient = np.divide(
        frozen_limited, frozen_difference,
        out=np.zeros_like(frozen_limited), where=frozen_difference != 0.0,
    )
    metrics["frozen_support"] = {
        "delta_energy": float(frozen_after - frozen_before),
        "delta_energy_over_dt": float((frozen_after - frozen_before) / dt),
        "gradient_dot_velocity": float(np.sum(frozen_gradient * frozen_rate)),
        "predicted_pair_dissipation": -float(np.sum(frozen_coefficient * frozen_difference**2)),
        "minimum_donor_factor": float(np.min(frozen_donor)),
    }
    return result, audit, metrics, (ecell, ei, ej), frozen_gradient


def main() -> int:
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("actual-state reproduction must run in a Slurm allocation")
    output = Path(sys.argv[1] if len(sys.argv) > 1 else "output/first-transition")
    output.mkdir(parents=True, exist_ok=True)
    config = cfg()
    eta0, seeds, orientations = voronoi_polycrystal(SHAPE, GRAINS, SEED, width=2.0, periodic=True)
    solver = MultiphaseFieldSolver(eta0, config, orientations=orientations)
    energies = [solver._anisotropic_energy()]
    reproduction = []
    captured = None
    for step in range(64):
        pre = solver.eta.copy()
        pre_active = solver.active_phases.copy()
        stats, _, _ = support_statistics(pre, pre_active)
        started = time.perf_counter()
        solver.step(compute_energy=False)
        elapsed = time.perf_counter() - started
        after_energy = solver._anisotropic_energy()
        delta = after_energy - energies[-1]
        energies.append(after_energy)
        record = {"step": step + 1, "energy_before": energies[-2], "energy_after": after_energy, "delta_energy": delta, "solver_seconds": elapsed, **stats}
        reproduction.append(record)
        (output / "reproduction-progress.json").write_text(json.dumps(reproduction, indent=2) + "\n")
        print(f"step={step + 1} E={after_energy:.17g} dE={delta:.17g} seconds={elapsed:.6g}", flush=True)
        tolerance = max(1e-10, abs(energies[-2]) * 1e-10)
        if delta > tolerance:
            captured = (step + 1, pre, solver.eta.copy(), pre_active, after_energy, delta, tolerance)
            break
    if captured is None:
        raise RuntimeError("no positive transition reproduced within 64 steps")

    step, pre, after, active, after_energy, observed_delta, tolerance = captured
    base_energy_support, base_update_support = support_masks(pre, active, True)
    trials = []
    primary = None
    for divisor in (1, 2, 4, 8, 16):
        result, audit, metrics, energy_edges, frozen_gradient = trial(
            pre, active, orientations, DT / divisor, config,
            (base_energy_support, base_update_support),
        )
        trials.append(metrics)
        if divisor == 1:
            primary = (result, audit, energy_edges, frozen_gradient)
    assert primary is not None
    result, audit, energy_edges, frozen_gradient = primary
    chemical, pair_rate, donor, local_count, cell, pi, pj, difference, raw, limited, coefficient, magnitude = audit
    executed_velocity = (result - pre) / DT
    before_components = anisotropic_energy_components(*energy_args(pre, active, orientations, config))
    after_components = anisotropic_energy_components(*energy_args(result, active, orientations, config))

    production_energy, production_gradient = anisotropic_energy_gradient(*energy_args(pre, active, orientations, config))
    oracle_energy, oracle_gradient = dense_streaming_oracle(
        pre, active, orientations, *law_args(config),
    )
    ecell, ei, ej = energy_edges
    frozen_energy = frozen_energy_gradient(pre, orientations, ecell, ei, ej, *law_args(config))[0]
    frozen_gdot = float(np.sum(frozen_gradient * executed_velocity))
    directional = []
    for factor in EPSILON_FACTORS:
        epsilon = DT * factor
        perturbed = pre + epsilon * executed_velocity
        production_perturbed = anisotropic_energy_gradient(*energy_args(perturbed, active, orientations, config))[0]
        oracle_perturbed = dense_streaming_oracle(
            perturbed, active, orientations, *law_args(config),
        )[0]
        frozen_perturbed = frozen_energy_gradient(
            perturbed, orientations, ecell, ei, ej, *law_args(config),
        )[0]
        directional.append({
            "epsilon": epsilon,
            "production_quotient": float((production_perturbed - production_energy) / epsilon),
            "dense_oracle_quotient": float((oracle_perturbed - oracle_energy) / epsilon),
            "frozen_support_quotient": float((frozen_perturbed - frozen_energy) / epsilon),
        })

    graph = graph_metrics(cell, pi, pj, difference, raw, limited, coefficient, local_count, donor)
    gdot = float(np.sum(production_gradient * executed_velocity))
    oracle_gdot = float(np.sum(oracle_gradient * executed_velocity))
    predicted = -float(np.sum(coefficient * difference**2))
    cell_gdot = np.sum(production_gradient * executed_velocity, axis=0)
    cell_dissipation = np.bincount(cell, weights=-coefficient * difference**2, minlength=SHAPE[0] * SHAPE[1]).reshape(SHAPE)
    pair_code = pi.astype(np.int64) * GRAINS + pj.astype(np.int64)
    pair_dissipation = np.bincount(pair_code, weights=-coefficient * difference**2, minlength=GRAINS * GRAINS).reshape(GRAINS, GRAINS)
    count_dissipation = np.bincount(local_count.ravel(), weights=cell_dissipation.ravel(), minlength=GRAINS + 1)
    limited_cell = np.bincount(cell, weights=(limited != raw).astype(np.int64), minlength=SHAPE[0] * SHAPE[1]).reshape(SHAPE)
    limited_mask = limited != raw
    active_count_summary = []
    for count in np.flatnonzero(np.bincount(local_count.ravel(), minlength=GRAINS + 1)):
        cells_at_count = local_count == count
        active_count_summary.append({
            "active_phase_count": int(count),
            "local_graph_degree": int(max(count - 1, 0)),
            "cell_count": int(np.count_nonzero(cells_at_count)),
            "gradient_dot_velocity": float(np.sum(cell_gdot[cells_at_count])),
            "predicted_pair_dissipation": float(np.sum(cell_dissipation[cells_at_count])),
            "donor_limited_edge_count": int(np.sum(limited_cell[cells_at_count])),
        })
    pair_values = pair_dissipation.ravel()
    strongest_pairs = np.argsort(np.abs(pair_values))[-20:][::-1]
    pair_summary = [
        {
            "phase_i": int(code // GRAINS),
            "phase_j": int(code % GRAINS),
            "predicted_pair_dissipation": float(pair_values[code]),
        }
        for code in strongest_pairs if pair_values[code] != 0.0
    ]
    cell_values = cell_gdot.ravel()
    largest_cells = np.argsort(cell_values)[-20:][::-1]
    cell_summary = [
        {
            "y": int(code // SHAPE[1]),
            "x": int(code % SHAPE[1]),
            "gradient_dot_velocity": float(cell_values[code]),
            "predicted_pair_dissipation": float(cell_dissipation.ravel()[code]),
            "active_phase_count": int(local_count.ravel()[code]),
            "local_graph_degree": int(max(local_count.ravel()[code] - 1, 0)),
            "donor_limited_edge_count": int(limited_cell.ravel()[code]),
        }
        for code in largest_cells
    ]

    relative_gradient_error = float(np.linalg.norm(production_gradient - oracle_gradient) / max(np.linalg.norm(oracle_gradient), np.finfo(float).tiny))
    directional_relative_error = abs(
        directional[-1]["production_quotient"] - gdot
    ) / max(1.0, abs(gdot))
    frozen_directional_relative_error = abs(
        directional[-1]["frozen_support_quotient"] - frozen_gdot
    ) / max(1.0, abs(frozen_gdot))
    if relative_gradient_error > 1e-10:
        classification = "SPARSE_SUPPORT_MISMATCH"
    elif gdot > 0.0 or predicted > 0.0 or graph["coefficient_minimum"] < 0.0:
        classification = "PAIR_LIMITER_NOT_SPD"
    elif directional_relative_error <= 1e-3:
        classification = "DERIVATIVE_MATCHES"
    elif frozen_directional_relative_error <= 1e-3:
        classification = "ACTIVE_BRANCH_NONDIFFERENTIABLE"
    elif (
        trials[0]["inclination_below_update_angle_tolerance"] > 0
        and abs(directional[-1]["production_quotient"] - gdot) > abs(gdot) * 1e-3
    ):
        classification = "LOW_GRADIENT_ORIENTATION_SINGULARITY"
    else:
        classification = "OTHER_IDENTIFIED_DEFECT"
    if classification == "DERIVATIVE_MATCHES" and trials[-1]["delta_energy"] < 0.0:
        root_cause = "MISSING_DENSE_GRAPH_TIMESTEP_FACTOR"
    elif trials[0]["frozen_support"]["delta_energy"] < 0.0 and trials[0]["delta_energy"] > 0.0:
        root_cause = "PAIR_SUPPORT_SWITCHING"
    elif classification == "SPARSE_SUPPORT_MISMATCH":
        root_cause = "SPARSE_SUPPORT_BOOKKEEPING"
    elif classification == "LOW_GRADIENT_ORIENTATION_SINGULARITY":
        root_cause = "LOW_GRADIENT_INCLINATION_EVALUATION"
    elif classification == "PAIR_LIMITER_NOT_SPD":
        root_cause = "PAIR_LIMITER_MOBILITY_OPERATOR"
    else:
        root_cause = "UNRESOLVED_FROM_DECLARED_DIAGNOSTICS"

    np.savez_compressed(
        output / "transition-state.npz", eta_before=pre, eta_after=after,
        active=active, energy_support=base_energy_support,
        update_support=base_update_support, orientations=orientations,
        seeds=seeds, chemical_potential=chemical,
        executed_velocity=executed_velocity, pair_rate=pair_rate,
        donor_factors=donor, local_active_count=local_count,
    )
    np.savez_compressed(
        output / "executed-pair-graph.npz", cell=cell, phase_i=pi, phase_j=pj,
        chemical_difference=difference, raw_exchange=raw,
        limited_exchange=limited, shared_coefficient=coefficient,
        pair_gradient_magnitude=magnitude,
    )
    np.savez_compressed(
        output / "audit-decomposition.npz", cell_gradient_dot_velocity=cell_gdot,
        cell_predicted_dissipation=cell_dissipation,
        phase_pair_predicted_dissipation=pair_dissipation,
        active_count_predicted_dissipation=count_dissipation,
        donor_limited_edges_per_cell=limited_cell,
        energy_potential_before=before_components[0], energy_gradient_before=before_components[1],
        energy_total_before=before_components[2], energy_potential_after=after_components[0],
        energy_gradient_after=after_components[1], energy_total_after=after_components[2],
        local_active_count=local_count,
    )
    report = {
        "schema": "anisotropic-first-positive-transition-audit-v1",
        "job_id": os.environ["SLURM_JOB_ID"],
        "configuration": {"shape": SHAPE, "grains": GRAINS, "seed": SEED, "accepted_dt": DT, "strength": "A2_STRONG", "anisotropic_energy": True, "anisotropic_mobility": False, "energy_normalization": ENERGY_NORMALIZATION, "mobility_normalization": MOBILITY_NORMALIZATION},
        "first_positive_transition": {"step": step, "energy_before": production_energy, "energy_after": after_energy, "delta_energy": observed_delta, "roundoff_tolerance": tolerance},
        "thermodynamic_audit": {
            "production_energy": production_energy,
            "dense_oracle_energy": oracle_energy,
            "frozen_support_energy": frozen_energy,
            "production_gradient_dot_executed_velocity": gdot,
            "dense_oracle_gradient_dot_executed_velocity": oracle_gdot,
            "frozen_support_gradient_dot_executed_velocity": frozen_gdot,
            "predicted_pair_dissipation": predicted,
            "production_dense_relative_gradient_error": relative_gradient_error,
            "production_directional_relative_error": directional_relative_error,
            "frozen_directional_relative_error": frozen_directional_relative_error,
            "executed_pair_graph": graph,
            "single_shared_nonnegative_coefficient": bool(
                graph["coefficient_minimum"] >= 0.0
                and graph["pair_representation_maximum_error"] <= 1e-14
            ),
            "decomposition": {
                "cell_largest_gradient_dot_velocity": cell_summary,
                "phase_pairs_largest_absolute_dissipation": pair_summary,
                "by_active_phase_count_and_graph_degree": active_count_summary,
                "donor_limiter": {
                    "unlimited_predicted_dissipation": float(np.sum((-coefficient * difference**2)[~limited_mask])),
                    "limited_predicted_dissipation": float(np.sum((-coefficient * difference**2)[limited_mask])),
                    "unlimited_edge_count": int(np.count_nonzero(~limited_mask)),
                    "limited_edge_count": int(np.count_nonzero(limited_mask)),
                },
            },
            "potential_energy_before": float(np.sum(before_components[0])),
            "gradient_energy_before": float(np.sum(before_components[1])),
            "potential_energy_after": float(np.sum(after_components[0])),
            "gradient_energy_after": float(np.sum(after_components[1])),
        },
        "directional_derivative": directional,
        "timestep_trials": trials,
        "support_proliferation": reproduction,
        "classification": classification,
        "root_cause": root_cause,
    }
    (output / "audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
