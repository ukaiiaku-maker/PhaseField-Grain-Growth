"""Native pairwise variational anisotropic phase-field operator.

The active-set count is held fixed while taking each discrete variation.  On
that smooth branch, ``anisotropic_energy_gradient`` is the exact derivative of
the energy it returns.  Pair exchange of the resulting chemical potentials is
antisymmetric, so the unconstrained update conserves the phase sum.
"""
from __future__ import annotations

import numpy as np
from numba import njit
from numpy.typing import NDArray

Array = NDArray[np.float64]


@njit(cache=True)
def _support_value_and_first(psi: float, power: int) -> tuple[float, float]:
    c = np.cos(psi)
    s = np.sin(psi)
    u = c ** power + s ** power
    value = u ** (1.0 / power)
    first = u ** (1.0 / power - 1.0) * (
        s ** (power - 1) * c - c ** (power - 1) * s
    )
    return value, first


@njit(cache=True)
def _pair_law(
    theta: float,
    phi_i: float,
    phi_j: float,
    gamma0: float,
    mobility0: float,
    g_min: float,
    inclination_weight: float,
    support_power: int,
    mobility_exponent: float,
    angular_scale: float,
    energy_normalization: float,
    mobility_normalization: float,
) -> tuple[float, float, float]:
    delta = abs((phi_i - phi_j + np.pi / 4.0) % (np.pi / 2.0) - np.pi / 4.0)
    misorientation = g_min + (1.0 - g_min) * np.sin(2.0 * delta) ** 2
    hi, dhi = _support_value_and_first(theta - phi_i, support_power)
    hj, dhj = _support_value_and_first(theta - phi_j, support_power)
    prefactor = gamma0 * energy_normalization * misorientation * angular_scale
    gamma = prefactor * (
        1.0 - inclination_weight
        + 0.5 * inclination_weight * (hi + hj)
    )
    gamma_theta = 0.5 * prefactor * inclination_weight * (dhi + dhj)
    mobility = (
        mobility0
        * mobility_normalization
        * (gamma / gamma0) ** (-mobility_exponent)
    )
    return gamma, gamma_theta, mobility


@njit(cache=True)
def _rotated_support_norm_and_gradient(
    px: float, py: float, phi: float, power: int,
) -> tuple[float, float, float]:
    """Return a rotated even p-norm and its Cartesian gradient."""
    c = np.cos(phi)
    s = np.sin(phi)
    a = c * px + s * py
    b = -s * px + c * py
    total = a ** power + b ** power
    if total == 0.0:
        return 0.0, 0.0, 0.0
    scale = total ** (1.0 / power - 1.0)
    da = scale * a ** (power - 1)
    db = scale * b ** (power - 1)
    return total ** (1.0 / power), c * da - s * db, s * da + c * db


@njit(cache=True)
def _pair_homogeneous_norm_and_gradient(
    px: float,
    py: float,
    phi_i: float,
    phi_j: float,
    gamma0: float,
    g_min: float,
    inclination_weight: float,
    support_power: int,
    angular_scale: float,
    energy_normalization: float,
) -> tuple[float, float, float]:
    """Return the one-homogeneous pair tension and Cahn--Hoffman vector."""
    delta = abs((phi_i - phi_j + np.pi / 4.0) % (np.pi / 2.0) - np.pi / 4.0)
    misorientation = g_min + (1.0 - g_min) * np.sin(2.0 * delta) ** 2
    prefactor = gamma0 * energy_normalization * misorientation * angular_scale
    radius = np.sqrt(px * px + py * py)
    radial_x = px / radius if radius > 0.0 else 0.0
    radial_y = py / radius if radius > 0.0 else 0.0
    hi, hix, hiy = _rotated_support_norm_and_gradient(px, py, phi_i, support_power)
    hj, hjx, hjy = _rotated_support_norm_and_gradient(px, py, phi_j, support_power)
    value = prefactor * (
        (1.0 - inclination_weight) * radius
        + 0.5 * inclination_weight * (hi + hj)
    )
    first_x = prefactor * (
        (1.0 - inclination_weight) * radial_x
        + 0.5 * inclination_weight * (hix + hjx)
    )
    first_y = prefactor * (
        (1.0 - inclination_weight) * radial_y
        + 0.5 * inclination_weight * (hiy + hjy)
    )
    return value, first_x, first_y


@njit(cache=True)
def anisotropic_energy_gradient(
    eta: Array,
    active: NDArray[np.bool_],
    orientations: Array,
    gamma0: float,
    mobility0: float,
    width: float,
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
) -> tuple[float, Array]:
    """Return one discrete pair energy and its exact nodal derivative."""
    phases, height, width_pixels = eta.shape
    derivative = np.zeros_like(eta)
    energy = 0.0
    gradient_scale = width * width / (np.pi * np.pi)
    density_scale = 4.0 / width
    tolerance = 0.0
    local = np.empty(phases, dtype=np.int64)

    for y in range(height):
        yp = (y + 1) % height if periodic else min(y + 1, height - 1)
        for x in range(width_pixels):
            xp = (x + 1) % width_pixels if periodic else min(x + 1, width_pixels - 1)
            count = 0
            for phase in range(phases):
                if active[phase] and (
                    eta[phase, y, x] > tolerance
                    or eta[phase, yp, x] > tolerance
                    or eta[phase, y, xp] > tolerance
                ):
                    local[count] = phase
                    count += 1
            if count < 2:
                continue
            # Each supported pair contributes the same fixed double-obstacle
            # gradient coefficient.  A local ``1 / count`` factor makes the
            # functional discontinuous when a third phase enters or leaves a
            # junction stencil: the coefficient of every surviving pair then
            # changes even though its fields did not.  Pair support itself is
            # safe to skip because an omitted pair has identically zero local
            # density; its coefficient must not depend on other phases.
            pair_gradient_scale = gradient_scale
            for left in range(count - 1):
                i = local[left]
                ui = eta[i, y, x]
                gix = (eta[i, y, xp] - ui) / dx
                giy = (eta[i, yp, x] - ui) / dx
                for right in range(left + 1, count):
                    j = local[right]
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
                        gradient = pair_gradient_scale * (
                            norm * norm - sum_norm * sum_norm
                        ) / (4.0 * gamma0)
                        common_x = pair_gradient_scale * sum_norm * sum_norm_x / (2.0 * gamma0)
                        common_y = pair_gradient_scale * sum_norm * sum_norm_y / (2.0 * gamma0)
                        anisotropic_x = pair_gradient_scale * norm * norm_x / (2.0 * gamma0)
                        anisotropic_y = pair_gradient_scale * norm * norm_y / (2.0 * gamma0)
                        flux_ix = density_scale * (anisotropic_x - common_x)
                        flux_iy = density_scale * (anisotropic_y - common_y)
                        flux_jx = density_scale * (-anisotropic_x - common_x)
                        flux_jy = density_scale * (-anisotropic_y - common_y)
                    else:
                        potential = gamma0 * ui * uj
                        gradient = -gamma0 * pair_gradient_scale * (
                            gix * gjx + giy * gjy
                        )
                        flux_ix = -density_scale * gamma0 * pair_gradient_scale * gjx
                        flux_iy = -density_scale * gamma0 * pair_gradient_scale * gjy
                        flux_jx = -density_scale * gamma0 * pair_gradient_scale * gix
                        flux_jy = -density_scale * gamma0 * pair_gradient_scale * giy
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


@njit(cache=True)
def anisotropic_pairwise_step(
    eta: Array,
    active: NDArray[np.bool_],
    orientations: Array,
    mobility_scale: Array,
    external: Array,
    use_external: bool,
    dt: float,
    gamma0: float,
    mobility0: float,
    width: float,
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
) -> tuple[Array, float, Array]:
    """Advance all pair forces with one work-conjugate mobility per pair."""
    energy, chemical = anisotropic_energy_gradient(
        eta, active, orientations, gamma0, mobility0, width, dx, periodic,
        g_min, inclination_weight, support_power, mobility_exponent,
        angular_scale, energy_normalization, mobility_normalization,
        anisotropic_energy,
    )
    phases, height, width_pixels = eta.shape
    rate = np.zeros_like(eta)
    local = np.empty(phases, dtype=np.int64)
    outgoing = np.empty(phases, dtype=np.float64)
    donor_scale = np.empty(phases, dtype=np.float64)
    kinetic_scale = np.pi * np.pi / (4.0 * width * dx * dx)
    tolerance = 0.0
    for y in range(height):
        ym = (y - 1) % height if periodic else max(y - 1, 0)
        yp = (y + 1) % height if periodic else min(y + 1, height - 1)
        for x in range(width_pixels):
            xm = (x - 1) % width_pixels if periodic else max(x - 1, 0)
            xp = (x + 1) % width_pixels if periodic else min(x + 1, width_pixels - 1)
            count = 0
            for phase in range(phases):
                outgoing[phase] = 0.0
                donor_scale[phase] = 1.0
                if active[phase] and (
                    eta[phase, y, x] > tolerance
                    or eta[phase, ym, x] > tolerance
                    or eta[phase, yp, x] > tolerance
                    or eta[phase, y, xm] > tolerance
                    or eta[phase, y, xp] > tolerance
                ):
                    local[count] = phase
                    count += 1
            for left in range(count - 1):
                i = local[left]
                for right in range(left + 1, count):
                    j = local[right]
                    px = (
                        eta[i, y, xp] - eta[i, y, xm]
                        - eta[j, y, xp] + eta[j, y, xm]
                    ) / (2.0 * dx)
                    py = (
                        eta[i, yp, x] - eta[i, ym, x]
                        - eta[j, yp, x] + eta[j, ym, x]
                    ) / (2.0 * dx)
                    theta = np.arctan2(py, px) if px * px + py * py > tolerance else 0.0
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
                    if exchange > 0.0:
                        outgoing[j] += exchange
                    else:
                        outgoing[i] -= exchange

            # Enforce the obstacle as a conservative pair-flux constraint.
            # Scaling every downhill exchange from an overdrawn donor by the
            # same nonnegative factor keeps each effective pair coefficient
            # symmetric and the nodal mobility graph positive semidefinite.
            for index in range(count):
                phase = local[index]
                demand = dt * outgoing[phase]
                if demand > eta[phase, y, x] and demand > 0.0:
                    donor_scale[phase] = eta[phase, y, x] / demand

            for left in range(count - 1):
                i = local[left]
                for right in range(left + 1, count):
                    j = local[right]
                    px = (
                        eta[i, y, xp] - eta[i, y, xm]
                        - eta[j, y, xp] + eta[j, y, xm]
                    ) / (2.0 * dx)
                    py = (
                        eta[i, yp, x] - eta[i, ym, x]
                        - eta[j, yp, x] + eta[j, ym, x]
                    ) / (2.0 * dx)
                    theta = np.arctan2(py, px) if px * px + py * py > tolerance else 0.0
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
                    exchange *= donor_scale[j] if exchange > 0.0 else donor_scale[i]
                    rate[i, y, x] += exchange
                    rate[j, y, x] -= exchange

    result = eta + dt * rate
    # The limiter makes the trial feasible analytically.  Remove only negative
    # roundoff and close the phase sum on the largest component; a larger
    # violation indicates an implementation error rather than a projection
    # event that may be silently normalized away.
    roundoff_tolerance = 1e-12
    zero_tolerance = 1e-14
    for y in range(height):
        for x in range(width_pixels):
            total = 0.0
            largest = 0
            for phase in range(phases):
                if result[phase, y, x] < -roundoff_tolerance:
                    raise FloatingPointError("pair-flux obstacle constraint failed")
                if abs(result[phase, y, x]) <= zero_tolerance:
                    result[phase, y, x] = 0.0
                if result[phase, y, x] > result[largest, y, x]:
                    largest = phase
                total += result[phase, y, x]
            result[largest, y, x] += 1.0 - total
            if result[largest, y, x] < -roundoff_tolerance:
                raise FloatingPointError("phase-sum roundoff correction failed")
    return result, energy, chemical


@njit(cache=True)
def anisotropic_pairwise_audit(
    eta: Array,
    active: NDArray[np.bool_],
    orientations: Array,
    mobility_scale: Array,
    dt: float,
    gamma0: float,
    mobility0: float,
    width: float,
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
) -> tuple[Array, Array, Array, Array, Array, Array, Array, Array, Array, Array, Array, Array]:
    """Replay one production pair step and expose its executed graph.

    This diagnostic uses the production support predicate, pair law, donor
    limiter, and loop ordering. Each row in the returned edge arrays is one
    undirected cell-local pair; ``limited`` is added to phase ``i`` and
    subtracted from phase ``j``.
    """
    _, chemical = anisotropic_energy_gradient(
        eta, active, orientations, gamma0, mobility0, width, dx, periodic,
        g_min, inclination_weight, support_power, mobility_exponent,
        angular_scale, energy_normalization, mobility_normalization,
        anisotropic_energy,
    )
    phases, height, width_pixels = eta.shape
    tolerance = 0.0
    local = np.empty(phases, dtype=np.int64)
    local_count = np.zeros((height, width_pixels), dtype=np.int32)
    edge_total = 0
    for y in range(height):
        ym = (y - 1) % height if periodic else max(y - 1, 0)
        yp = (y + 1) % height if periodic else min(y + 1, height - 1)
        for x in range(width_pixels):
            xm = (x - 1) % width_pixels if periodic else max(x - 1, 0)
            xp = (x + 1) % width_pixels if periodic else min(x + 1, width_pixels - 1)
            count = 0
            for phase in range(phases):
                if active[phase] and (
                    eta[phase, y, x] > tolerance
                    or eta[phase, ym, x] > tolerance
                    or eta[phase, yp, x] > tolerance
                    or eta[phase, y, xm] > tolerance
                    or eta[phase, y, xp] > tolerance
                ):
                    count += 1
            local_count[y, x] = count
            edge_total += count * (count - 1) // 2

    cell = np.empty(edge_total, dtype=np.int32)
    phase_i = np.empty(edge_total, dtype=np.int16)
    phase_j = np.empty(edge_total, dtype=np.int16)
    mu_difference = np.empty(edge_total, dtype=np.float64)
    raw = np.empty(edge_total, dtype=np.float64)
    limited = np.empty(edge_total, dtype=np.float64)
    coefficient = np.empty(edge_total, dtype=np.float64)
    gradient_magnitude = np.empty(edge_total, dtype=np.float64)
    donor_factors = np.ones_like(eta)
    rate = np.zeros_like(eta)
    outgoing = np.empty(phases, dtype=np.float64)
    donor_scale = np.empty(phases, dtype=np.float64)
    kinetic_scale = np.pi * np.pi / (4.0 * width * dx * dx)
    edge = 0
    for y in range(height):
        ym = (y - 1) % height if periodic else max(y - 1, 0)
        yp = (y + 1) % height if periodic else min(y + 1, height - 1)
        for x in range(width_pixels):
            xm = (x - 1) % width_pixels if periodic else max(x - 1, 0)
            xp = (x + 1) % width_pixels if periodic else min(x + 1, width_pixels - 1)
            count = 0
            for phase in range(phases):
                outgoing[phase] = 0.0
                donor_scale[phase] = 1.0
                if active[phase] and (
                    eta[phase, y, x] > tolerance
                    or eta[phase, ym, x] > tolerance
                    or eta[phase, yp, x] > tolerance
                    or eta[phase, y, xm] > tolerance
                    or eta[phase, y, xp] > tolerance
                ):
                    local[count] = phase
                    count += 1
            begin = edge
            for left in range(count - 1):
                i = local[left]
                for right in range(left + 1, count):
                    j = local[right]
                    px = (
                        eta[i, y, xp] - eta[i, y, xm]
                        - eta[j, y, xp] + eta[j, y, xm]
                    ) / (2.0 * dx)
                    py = (
                        eta[i, yp, x] - eta[i, ym, x]
                        - eta[j, yp, x] + eta[j, ym, x]
                    ) / (2.0 * dx)
                    magnitude = np.sqrt(px * px + py * py)
                    theta = np.arctan2(py, px) if magnitude * magnitude > tolerance else 0.0
                    _, _, law_mobility = _pair_law(
                        theta, orientations[i], orientations[j], gamma0, mobility0,
                        g_min, inclination_weight, support_power, mobility_exponent,
                        angular_scale, energy_normalization, mobility_normalization,
                    )
                    pair_mobility = law_mobility if anisotropic_mobility else mobility0
                    difference = chemical[j, y, x] - chemical[i, y, x]
                    exchange = pair_mobility * kinetic_scale * difference * mobility_scale[y, x]
                    cell[edge] = y * width_pixels + x
                    phase_i[edge] = i
                    phase_j[edge] = j
                    mu_difference[edge] = difference
                    raw[edge] = exchange
                    gradient_magnitude[edge] = magnitude
                    if exchange > 0.0:
                        outgoing[j] += exchange
                    else:
                        outgoing[i] -= exchange
                    edge += 1
            for index in range(count):
                phase = local[index]
                demand = dt * outgoing[phase]
                if demand > eta[phase, y, x] and demand > 0.0:
                    donor_scale[phase] = eta[phase, y, x] / demand
                donor_factors[phase, y, x] = donor_scale[phase]
            for entry in range(begin, edge):
                i = int(phase_i[entry])
                j = int(phase_j[entry])
                exchange = raw[entry]
                scale = donor_scale[j] if exchange > 0.0 else donor_scale[i]
                executed = exchange * scale
                limited[entry] = executed
                coefficient[entry] = (
                    executed / mu_difference[entry]
                    if mu_difference[entry] != 0.0 else 0.0
                )
                rate[i, y, x] += executed
                rate[j, y, x] -= executed
    return (
        chemical, rate, donor_factors, local_count, cell, phase_i, phase_j,
        mu_difference, raw, limited, coefficient, gradient_magnitude,
    )


@njit(cache=True)
def anisotropic_energy_components(
    eta: Array,
    active: NDArray[np.bool_],
    orientations: Array,
    gamma0: float,
    mobility0: float,
    width: float,
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
) -> tuple[Array, Array, Array]:
    """Return cellwise potential, gradient, and total implemented energies."""
    phases, height, width_pixels = eta.shape
    potential = np.zeros((height, width_pixels), dtype=np.float64)
    gradient = np.zeros((height, width_pixels), dtype=np.float64)
    gradient_scale = width * width / (np.pi * np.pi)
    density_scale = 4.0 / width
    tolerance = 1e-14
    local = np.empty(phases, dtype=np.int64)
    for y in range(height):
        yp = (y + 1) % height if periodic else min(y + 1, height - 1)
        for x in range(width_pixels):
            xp = (x + 1) % width_pixels if periodic else min(x + 1, width_pixels - 1)
            count = 0
            for phase in range(phases):
                if active[phase] and (
                    eta[phase, y, x] > tolerance
                    or eta[phase, yp, x] > tolerance
                    or eta[phase, y, xp] > tolerance
                ):
                    local[count] = phase
                    count += 1
            for left in range(count - 1):
                i = local[left]
                ui = eta[i, y, x]
                gix = (eta[i, y, xp] - ui) / dx
                giy = (eta[i, yp, x] - ui) / dx
                for right in range(left + 1, count):
                    j = local[right]
                    uj = eta[j, y, x]
                    gjx = (eta[j, y, xp] - uj) / dx
                    gjy = (eta[j, yp, x] - uj) / dx
                    px = gix - gjx
                    py = giy - gjy
                    if anisotropic_energy:
                        norm, _, _ = _pair_homogeneous_norm_and_gradient(
                            px, py, orientations[i], orientations[j], gamma0,
                            g_min, inclination_weight, support_power, angular_scale,
                            energy_normalization,
                        )
                        sx = gix + gjx
                        sy = giy + gjy
                        sum_norm, _, _ = _pair_homogeneous_norm_and_gradient(
                            sx, sy, orientations[i], orientations[j], gamma0,
                            g_min, inclination_weight, support_power, angular_scale,
                            energy_normalization,
                        )
                        gradient_term = gradient_scale * (
                            norm * norm - sum_norm * sum_norm
                        ) / (4.0 * gamma0)
                    else:
                        gradient_term = -gamma0 * gradient_scale * (
                            gix * gjx + giy * gjy
                        )
                    potential[y, x] += density_scale * gamma0 * ui * uj * dx * dx
                    gradient[y, x] += density_scale * gradient_term * dx * dx
    return potential, gradient, potential + gradient
