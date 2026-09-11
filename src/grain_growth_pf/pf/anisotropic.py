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
            if count < 2:
                continue
            pair_gradient_scale = gradient_scale / count
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
                    magnitude = np.sqrt(px * px + py * py)
                    theta = np.arctan2(py, px) if magnitude > tolerance else 0.0
                    law_gamma, law_first, _ = _pair_law(
                        theta, orientations[i], orientations[j], gamma0, mobility0,
                        g_min, inclination_weight, support_power, mobility_exponent,
                        angular_scale, energy_normalization, mobility_normalization,
                    )
                    gamma = law_gamma if anisotropic_energy else gamma0
                    gamma_first = law_first if anisotropic_energy else 0.0
                    cross = gix * gjx + giy * gjy
                    base = ui * uj - pair_gradient_scale * cross
                    energy += density_scale * gamma * base * dx * dx

                    dgamma_x = 0.0
                    dgamma_y = 0.0
                    if magnitude > tolerance:
                        dgamma_x = -gamma_first * py / (magnitude * magnitude)
                        dgamma_y = gamma_first * px / (magnitude * magnitude)
                    flux_ix = density_scale * (
                        -gamma * pair_gradient_scale * gjx + base * dgamma_x
                    )
                    flux_iy = density_scale * (
                        -gamma * pair_gradient_scale * gjy + base * dgamma_y
                    )
                    flux_jx = density_scale * (
                        -gamma * pair_gradient_scale * gix - base * dgamma_x
                    )
                    flux_jy = density_scale * (
                        -gamma * pair_gradient_scale * giy - base * dgamma_y
                    )
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
    kinetic_scale = np.pi * np.pi / (4.0 * width * dx * dx)
    tolerance = 1e-14
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
                    rate[i, y, x] += exchange
                    rate[j, y, x] -= exchange

    result = eta + dt * rate
    for y in range(height):
        for x in range(width_pixels):
            total = 0.0
            for phase in range(phases):
                if result[phase, y, x] < 0.0:
                    result[phase, y, x] = 0.0
                elif result[phase, y, x] > 1.0:
                    result[phase, y, x] = 1.0
                total += result[phase, y, x]
            if total > 0.0:
                for phase in range(phases):
                    result[phase, y, x] /= total
    return result, energy, chemical
