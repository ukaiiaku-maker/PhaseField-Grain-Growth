from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def circular_grain(shape: tuple[int, int], radius: float, width: float,
                   center: tuple[float, float] | None = None) -> NDArray[np.float64]:
    """Two-order-parameter circle with a compact double-obstacle interface."""
    center = center or ((shape[0] - 1) / 2, (shape[1] - 1) / 2)
    y, x = np.indices(shape, dtype=float)
    r = np.hypot(y - center[0], x - center[1])
    argument = np.pi * (r - radius) / width
    inside = np.where(
        argument <= -np.pi / 2, 1.0,
        np.where(argument >= np.pi / 2, 0.0, 0.5 * (1.0 - np.sin(argument))),
    )
    return np.stack((1.0 - inside, inside))


def planar_interface(shape: tuple[int, int], width: float, angle: float = 0.0,
                     offset: float = 0.0) -> NDArray[np.float64]:
    y, x = np.indices(shape, dtype=float)
    # Put an axis-aligned interface center on a grid point.  For the compact
    # obstacle profile this also puts its endpoints on grid points when the
    # configured width is an even number of cells, giving exact discrete
    # stationarity rather than a half-cell truncation error.
    yc, xc = np.asarray(shape, dtype=float) // 2
    signed = (x - xc) * np.cos(angle) + (y - yc) * np.sin(angle) - offset
    argument = np.pi * signed / width
    phase = np.where(
        argument <= -np.pi / 2, 0.0,
        np.where(argument >= np.pi / 2, 1.0, 0.5 * (1.0 + np.sin(argument))),
    )
    return np.stack((1.0 - phase, phase))


def voronoi_polycrystal(shape: tuple[int, int], n_grains: int, seed: int,
                        width: float = 2.0, periodic: bool = True) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    rng = np.random.default_rng(seed)
    seeds = rng.uniform([0, 0], shape, size=(n_grains, 2))
    orientations = rng.uniform(0.0, np.pi, n_grains)
    y, x = np.indices(shape, dtype=float)
    distances = []
    for sy, sx in seeds:
        dy, dx = np.abs(y - sy), np.abs(x - sx)
        if periodic:
            dy, dx = np.minimum(dy, shape[0] - dy), np.minimum(dx, shape[1] - dx)
        distances.append(dx * dx + dy * dy)
    distance = np.stack(distances)
    # Smooth hard Voronoi cells over a compact band.  The previous global
    # softmax left exponentially small tails for every grain at every pixel;
    # those tails make the local phase count ill-defined in a pairwise MPF.
    from scipy.ndimage import gaussian_filter

    labels = np.argmin(np.stack(distances), axis=0)
    eta = np.eye(n_grains, dtype=float)[labels].transpose(2, 0, 1)
    sigma = max(width / 2.0, 0.25)
    spatial_mode = "wrap" if periodic else "nearest"
    eta = gaussian_filter(
        eta, sigma=(0.0, sigma, sigma),
        mode=("nearest", spatial_mode, spatial_mode), truncate=2.0,
    )
    eta[eta < 1e-14] = 0.0
    eta /= eta.sum(axis=0, keepdims=True)
    return eta, seeds, orientations


def equivalent_radius(phase: NDArray[np.float64], dx: float = 1.0) -> float:
    return float(np.sqrt(phase.sum() * dx**2 / np.pi))
