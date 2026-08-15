from __future__ import annotations

import numpy as np
from numba import njit
from numpy.typing import NDArray


Array = NDArray[np.float64]


@njit(cache=True)
def pairwise_obstacle_step(
    eta: Array,
    active: NDArray[np.bool_],
    mobility_scale: Array,
    external: Array,
    use_external: bool,
    dt: float,
    mobility: float,
    gamma: float,
    width: float,
    dx: float,
    periodic: bool,
) -> Array:
    """Advance the Qiu pairwise equation without dense phase temporaries.

    The loops deliberately retain the same nine-point Laplacian and local
    one-cardinal-cell support as the vectorized reference formulation.  A
    compiled dense scan is cheaper than allocating several ``P x H x W``
    work arrays and immediately skips the compact bulk of every phase.
    """
    phases, height, width_pixels = eta.shape
    phase_sum = np.zeros((height, width_pixels), dtype=np.float64)
    lap_sum = np.zeros((height, width_pixels), dtype=np.float64)
    external_sum = np.zeros((height, width_pixels), dtype=np.float64)
    count = np.zeros((height, width_pixels), dtype=np.int32)
    inverse_lap_scale = 1.0 / (6.0 * dx * dx)

    for phase in range(phases):
        if not active[phase]:
            continue
        for y in range(height):
            ym = (y - 1) % height if periodic else max(y - 1, 0)
            yp = (y + 1) % height if periodic else min(y + 1, height - 1)
            for x in range(width_pixels):
                xm = (x - 1) % width_pixels if periodic else max(x - 1, 0)
                xp = (x + 1) % width_pixels if periodic else min(x + 1, width_pixels - 1)
                if not (
                    eta[phase, y, x] > 1e-14
                    or eta[phase, ym, x] > 1e-14
                    or eta[phase, yp, x] > 1e-14
                    or eta[phase, y, xm] > 1e-14
                    or eta[phase, y, xp] > 1e-14
                ):
                    continue
                value = eta[phase, y, x]
                lap = (
                    4.0
                    * (
                        eta[phase, ym, x]
                        + eta[phase, yp, x]
                        + eta[phase, y, xm]
                        + eta[phase, y, xp]
                    )
                    + eta[phase, ym, xm]
                    + eta[phase, ym, xp]
                    + eta[phase, yp, xm]
                    + eta[phase, yp, xp]
                    - 20.0 * value
                ) * inverse_lap_scale
                phase_sum[y, x] += value
                lap_sum[y, x] += lap
                if use_external:
                    external_sum[y, x] += external[phase, y, x]
                count[y, x] += 1

    obstacle = 2.0 * np.sin(np.pi * dx / (2.0 * width)) ** 2 / (dx * dx)
    result = eta.copy()
    trial_sum = np.zeros((height, width_pixels), dtype=np.float64)
    for phase in range(phases):
        if not active[phase]:
            continue
        for y in range(height):
            ym = (y - 1) % height if periodic else max(y - 1, 0)
            yp = (y + 1) % height if periodic else min(y + 1, height - 1)
            for x in range(width_pixels):
                xm = (x - 1) % width_pixels if periodic else max(x - 1, 0)
                xp = (x + 1) % width_pixels if periodic else min(x + 1, width_pixels - 1)
                if not (
                    eta[phase, y, x] > 1e-14
                    or eta[phase, ym, x] > 1e-14
                    or eta[phase, yp, x] > 1e-14
                    or eta[phase, y, xm] > 1e-14
                    or eta[phase, y, xp] > 1e-14
                ):
                    continue
                value = eta[phase, y, x]
                lap = (
                    4.0
                    * (
                        eta[phase, ym, x]
                        + eta[phase, yp, x]
                        + eta[phase, y, xm]
                        + eta[phase, y, xp]
                    )
                    + eta[phase, ym, xm]
                    + eta[phase, ym, xp]
                    + eta[phase, yp, xm]
                    + eta[phase, yp, xp]
                    - 20.0 * value
                ) * inverse_lap_scale
                rate = mobility * gamma * (
                    lap * phase_sum[y, x]
                    - value * lap_sum[y, x]
                    + obstacle * (count[y, x] * value - phase_sum[y, x])
                )
                if use_external:
                    rate += mobility * (
                        external[phase, y, x] - external_sum[y, x] / count[y, x]
                    )
                trial = value + dt * rate * mobility_scale[y, x]
                if trial < 0.0:
                    trial = 0.0
                elif trial > 1.0:
                    trial = 1.0
                result[phase, y, x] = trial
                trial_sum[y, x] += trial

    for phase in range(phases):
        if not active[phase]:
            continue
        for y in range(height):
            for x in range(width_pixels):
                if trial_sum[y, x] > 0.0:
                    result[phase, y, x] /= trial_sum[y, x]
    return result


@njit(cache=True)
def pairwise_free_energy(
    eta: Array,
    gamma: float,
    interface_width: float,
    dx: float,
    periodic: bool,
) -> float:
    """Evaluate the pairwise double-obstacle energy without work arrays."""
    phases, height, width_pixels = eta.shape
    density_sum = 0.0
    gradient_scale = interface_width * interface_width / (np.pi * np.pi)
    for y in range(height):
        yp = (y + 1) % height if periodic else min(y + 1, height - 1)
        for x in range(width_pixels):
            xp = (x + 1) % width_pixels if periodic else min(x + 1, width_pixels - 1)
            phase_sum = 0.0
            square_sum = 0.0
            gx_sum = 0.0
            gy_sum = 0.0
            gradient_square_sum = 0.0
            for phase in range(phases):
                value = eta[phase, y, x]
                gx = (eta[phase, yp, x] - value) / dx if periodic or y < height - 1 else 0.0
                gy = (eta[phase, y, xp] - value) / dx if periodic or x < width_pixels - 1 else 0.0
                phase_sum += value
                square_sum += value * value
                gx_sum += gx
                gy_sum += gy
                gradient_square_sum += gx * gx + gy * gy
            pair_potential = 0.5 * (phase_sum * phase_sum - square_sum)
            pair_gradient = 0.5 * (
                gx_sum * gx_sum + gy_sum * gy_sum - gradient_square_sum
            )
            density_sum += pair_potential - gradient_scale * pair_gradient
    return 4.0 * gamma / interface_width * density_sum * dx * dx
