"""Compact-support constrained update for the anisotropic PF energy.

The step first evaluates the exact discrete energy derivative, constructs a
cell-local candidate graph from exact support in the von Neumann stencil, and
forms the usual antisymmetric pair-mobility descent direction.  The trial is
then projected onto the simplex on that candidate graph.  Projection is the
solution of the local double-obstacle variational inequality and therefore
produces exact zeros.  A deterministic backtracking search accepts only an
unforced energy-decreasing step.

The candidate graph is a pure function of the checkpointed phase field.  A
phase absent from the current cell and its stencil cannot enter the cell.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit
from numpy.typing import NDArray

from .anisotropic import _pair_homogeneous_norm_and_gradient, _pair_law

Array = NDArray[np.float64]


@dataclass(frozen=True)
class CompactSupportDiagnostics:
    active_mean: float
    active_p95: float
    active_max: int
    candidate_mean: float
    candidate_p95: float
    candidate_max: int
    active_pair_instances: int
    candidate_pair_instances: int
    kkt_residual: float
    entries: int
    retirements: int
    line_search_scale: float


def exact_candidate_graph(
    eta: Array, periodic: bool,
) -> tuple[NDArray[np.int64], NDArray[np.int32]]:
    """Return padded cell candidates from exact current/stencil support.

    Phase indices are sorted, making the graph deterministic across restart.
    The vectorized support enumeration avoids a Python-level global-phase scan
    at every cell; all subsequent pair work is proportional to local support.
    """
    support = np.moveaxis(np.asarray(eta) > 0.0, 0, -1)
    candidate = support.copy()
    if periodic:
        candidate |= np.roll(support, 1, axis=0)
        candidate |= np.roll(support, -1, axis=0)
        candidate |= np.roll(support, 1, axis=1)
        candidate |= np.roll(support, -1, axis=1)
    else:
        candidate[1:] |= support[:-1]
        candidate[:-1] |= support[1:]
        candidate[:, 1:] |= support[:, :-1]
        candidate[:, :-1] |= support[:, 1:]
    counts = candidate.sum(axis=2, dtype=np.int32)
    maximum = max(1, int(np.max(counts)))
    graph = np.full((*counts.shape, maximum), -1, dtype=np.int64)
    rows, columns, phases = np.nonzero(candidate)
    cursor = np.zeros_like(counts)
    for y, x, phase in zip(rows, columns, phases, strict=True):
        slot = cursor[y, x]
        graph[y, x, slot] = phase
        cursor[y, x] += 1
    return graph, counts


@njit(cache=True)
def _candidate_energy_gradient(
    eta: Array,
    graph: NDArray[np.int64],
    counts: NDArray[np.int32],
    orientations: Array,
    gamma0: float,
    interface_width: float,
    dx: float,
    periodic: bool,
    g_min: float,
    inclination_weight: float,
    support_power: int,
    angular_scale: float,
    energy_normalization: float,
    anisotropic_energy: bool,
) -> tuple[float, Array]:
    """Evaluate energy and its derivative on one frozen local graph."""
    _, height, width_pixels = eta.shape
    derivative = np.zeros_like(eta)
    energy = 0.0
    gradient_scale = interface_width * interface_width / (np.pi * np.pi)
    density_scale = 4.0 / interface_width
    for y in range(height):
        yp = (y + 1) % height if periodic else min(y + 1, height - 1)
        for x in range(width_pixels):
            xp = (x + 1) % width_pixels if periodic else min(x + 1, width_pixels - 1)
            count = counts[y, x]
            for left in range(count - 1):
                i = graph[y, x, left]
                ui = eta[i, y, x]
                gix = (eta[i, y, xp] - ui) / dx
                giy = (eta[i, yp, x] - ui) / dx
                for right in range(left + 1, count):
                    j = graph[y, x, right]
                    uj = eta[j, y, x]
                    gjx = (eta[j, y, xp] - uj) / dx
                    gjy = (eta[j, yp, x] - uj) / dx
                    px = gix - gjx
                    py = giy - gjy
                    if anisotropic_energy:
                        norm, norm_x, norm_y = _pair_homogeneous_norm_and_gradient(
                            px, py, orientations[i], orientations[j], gamma0,
                            g_min, inclination_weight, support_power, angular_scale,
                            energy_normalization,
                        )
                        sx = gix + gjx
                        sy = giy + gjy
                        sum_norm, sum_norm_x, sum_norm_y = _pair_homogeneous_norm_and_gradient(
                            sx, sy, orientations[i], orientations[j], gamma0,
                            g_min, inclination_weight, support_power, angular_scale,
                            energy_normalization,
                        )
                        potential = gamma0 * ui * uj
                        gradient = gradient_scale * (norm * norm - sum_norm * sum_norm) / (4.0 * gamma0)
                        common_x = gradient_scale * sum_norm * sum_norm_x / (2.0 * gamma0)
                        common_y = gradient_scale * sum_norm * sum_norm_y / (2.0 * gamma0)
                        anisotropic_x = gradient_scale * norm * norm_x / (2.0 * gamma0)
                        anisotropic_y = gradient_scale * norm * norm_y / (2.0 * gamma0)
                        flux_ix = density_scale * (anisotropic_x - common_x)
                        flux_iy = density_scale * (anisotropic_y - common_y)
                        flux_jx = density_scale * (-anisotropic_x - common_x)
                        flux_jy = density_scale * (-anisotropic_y - common_y)
                    else:
                        potential = gamma0 * ui * uj
                        gradient = -gamma0 * gradient_scale * (gix * gjx + giy * gjy)
                        flux_ix = -density_scale * gamma0 * gradient_scale * gjx
                        flux_iy = -density_scale * gamma0 * gradient_scale * gjy
                        flux_jx = -density_scale * gamma0 * gradient_scale * gix
                        flux_jy = -density_scale * gamma0 * gradient_scale * giy
                    energy += density_scale * (potential + gradient) * dx * dx
                    derivative[i, y, x] += density_scale * gamma0 * uj * dx * dx
                    derivative[j, y, x] += density_scale * gamma0 * ui * dx * dx
                    derivative[i, y, x] -= (flux_ix + flux_iy) * dx
                    derivative[i, y, xp] += flux_ix * dx
                    derivative[i, yp, x] += flux_iy * dx
                    derivative[j, y, x] -= (flux_jx + flux_jy) * dx
                    derivative[j, y, xp] += flux_jx * dx
                    derivative[j, yp, x] += flux_jy * dx
    return energy, derivative


def compact_support_energy(
    eta: Array, orientations: Array, gamma0: float, interface_width: float,
    dx: float, periodic: bool, g_min: float, inclination_weight: float,
    support_power: int, angular_scale: float, energy_normalization: float,
    anisotropic_energy: bool,
) -> float:
    """Return the canonical energy on the deterministic local graph."""
    graph, counts = exact_candidate_graph(eta, periodic)
    return float(_candidate_energy_gradient(
        eta, graph, counts, orientations, gamma0, interface_width, dx,
        periodic, g_min, inclination_weight, support_power, angular_scale,
        energy_normalization, anisotropic_energy,
    )[0])


@njit(cache=True)
def _candidate_pair_rate(
    eta: Array,
    chemical: Array,
    graph: NDArray[np.int64],
    counts: NDArray[np.int32],
    orientations: Array,
    mobility_scale: Array,
    external: Array,
    use_external: bool,
    gamma0: float,
    mobility0: float,
    interface_width: float,
    dx: float,
    periodic: bool,
    g_min: float,
    inclination_weight: float,
    support_power: int,
    mobility_exponent: float,
    angular_scale: float,
    energy_normalization: float,
    mobility_normalization: float,
    anisotropic_mobility: bool,
) -> Array:
    phases, height, width = eta.shape
    rate = np.zeros_like(eta)
    kinetic_scale = np.pi * np.pi / (4.0 * interface_width * dx * dx)
    for y in range(height):
        ym = (y - 1) % height if periodic else max(y - 1, 0)
        yp = (y + 1) % height if periodic else min(y + 1, height - 1)
        for x in range(width):
            xm = (x - 1) % width if periodic else max(x - 1, 0)
            xp = (x + 1) % width if periodic else min(x + 1, width - 1)
            count = counts[y, x]
            for left in range(count - 1):
                i = graph[y, x, left]
                for right in range(left + 1, count):
                    j = graph[y, x, right]
                    px = (
                        eta[i, y, xp] - eta[i, y, xm]
                        - eta[j, y, xp] + eta[j, y, xm]
                    ) / (2.0 * dx)
                    py = (
                        eta[i, yp, x] - eta[i, ym, x]
                        - eta[j, yp, x] + eta[j, ym, x]
                    ) / (2.0 * dx)
                    theta = np.arctan2(py, px) if px * px + py * py > 0.0 else 0.0
                    _, _, law_mobility = _pair_law(
                        theta, orientations[i], orientations[j], gamma0, mobility0,
                        g_min, inclination_weight, support_power, mobility_exponent,
                        angular_scale, energy_normalization, mobility_normalization,
                    )
                    pair_mobility = law_mobility if anisotropic_mobility else mobility0
                    drive = kinetic_scale * (chemical[j, y, x] - chemical[i, y, x])
                    if use_external:
                        drive += (external[i, y, x] - external[j, y, x]) / count
                    exchange = pair_mobility * drive * mobility_scale[y, x]
                    rate[i, y, x] += exchange
                    rate[j, y, x] -= exchange
    return rate


@njit(cache=True)
def _project_candidate_simplex(
    eta: Array,
    rate: Array,
    graph: NDArray[np.int64],
    counts: NDArray[np.int32],
    step_size: float,
) -> tuple[Array, float]:
    phases, height, width = eta.shape
    result = np.zeros_like(eta)
    maximum = graph.shape[2]
    values = np.empty(maximum, dtype=np.float64)
    ordered = np.empty(maximum, dtype=np.float64)
    max_residual = 0.0
    for y in range(height):
        for x in range(width):
            count = counts[y, x]
            for slot in range(count):
                phase = graph[y, x, slot]
                values[slot] = eta[phase, y, x] + step_size * rate[phase, y, x]
                ordered[slot] = values[slot]
            ordered[:count].sort()
            cumulative = 0.0
            rho = 0
            theta = 0.0
            for reverse in range(count):
                value = ordered[count - reverse - 1]
                cumulative += value
                trial_theta = (cumulative - 1.0) / (reverse + 1)
                if value - trial_theta > 0.0:
                    rho = reverse + 1
                    theta = trial_theta
            if rho == 0:
                raise FloatingPointError("empty simplex active set")
            total = 0.0
            largest_phase = graph[y, x, 0]
            for slot in range(count):
                phase = graph[y, x, slot]
                projected = max(values[slot] - theta, 0.0)
                result[phase, y, x] = projected
                total += projected
                if projected > result[largest_phase, y, x]:
                    largest_phase = phase
                multiplier = projected - values[slot] + theta
                max_residual = max(
                    max_residual,
                    max(-multiplier, 0.0),
                    abs(projected * multiplier),
                )
            result[largest_phase, y, x] += 1.0 - total
            max_residual = max(max_residual, abs(total - 1.0))
    return result, max_residual


def compact_support_step(
    eta: Array,
    active: NDArray[np.bool_],
    orientations: Array,
    mobility_scale: Array,
    external: Array,
    use_external: bool,
    dt: float,
    gamma0: float,
    mobility0: float,
    interface_width: float,
    dx: float,
    periodic: bool,
    g_min: float,
    inclination_weight: float,
    support_power: int,
    mobility_exponent: float,
    angular_scale: float,
    energy_normalization: float,
    mobility_normalization: float,
    anisotropic_energy: bool,
    anisotropic_mobility: bool,
    kkt_tolerance: float,
) -> tuple[Array, float, Array, CompactSupportDiagnostics]:
    """Take one exact-support obstacle step with deterministic backtracking."""
    if kkt_tolerance not in (1e-8, 1e-10, 1e-12):
        raise ValueError("compact-support KKT tolerance must be 1e-8, 1e-10, or 1e-12")
    graph, counts = exact_candidate_graph(eta, periodic)
    energy, chemical = _candidate_energy_gradient(
        eta, graph, counts, orientations, gamma0, interface_width, dx,
        periodic, g_min, inclination_weight, support_power, angular_scale,
        energy_normalization, anisotropic_energy,
    )
    rate = _candidate_pair_rate(
        eta, chemical, graph, counts, orientations, mobility_scale, external,
        use_external, gamma0, mobility0, interface_width, dx, periodic, g_min,
        inclination_weight, support_power, mobility_exponent, angular_scale,
        energy_normalization, mobility_normalization, anisotropic_mobility,
    )
    scale = 1.0
    accepted = None
    residual = float("inf")
    for _ in range(40):
        trial, residual = _project_candidate_simplex(
            eta, rate, graph, counts, dt * scale
        )
        if use_external:
            accepted = trial
            break
        trial_energy = _candidate_energy_gradient(
            trial, graph, counts, orientations, gamma0, interface_width, dx,
            periodic, g_min, inclination_weight, support_power, angular_scale,
            energy_normalization, anisotropic_energy,
        )[0]
        energy_slack = 64.0 * np.finfo(float).eps * max(1.0, abs(energy))
        if trial_energy <= energy + energy_slack:
            accepted = trial
            break
        scale *= 0.5
    if accepted is None:
        raise FloatingPointError("compact-support energy line search failed")

    old_support = eta > 0.0
    new_support = accepted > 0.0
    active_counts = new_support.sum(axis=0)
    active_pairs = active_counts * (active_counts - 1) // 2
    candidate_pairs = counts * (counts - 1) // 2
    diagnostics = CompactSupportDiagnostics(
        active_mean=float(np.mean(active_counts)),
        active_p95=float(np.percentile(active_counts, 95)),
        active_max=int(np.max(active_counts)),
        candidate_mean=float(np.mean(counts)),
        candidate_p95=float(np.percentile(counts, 95)),
        candidate_max=int(np.max(counts)),
        active_pair_instances=int(np.sum(active_pairs)),
        candidate_pair_instances=int(np.sum(candidate_pairs)),
        kkt_residual=float(residual),
        entries=int(np.count_nonzero(new_support & ~old_support)),
        retirements=int(np.count_nonzero(old_support & ~new_support)),
        line_search_scale=scale,
    )
    return accepted, float(energy), chemical, diagnostics
