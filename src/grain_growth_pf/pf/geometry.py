from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def circular_grain(shape: tuple[int, int], radius: float, width: float,
                   center: tuple[float, float] | None = None) -> NDArray[np.float64]:
    """Two-order-parameter circular grain with a tanh diffuse interface."""
    center = center or ((shape[0] - 1) / 2, (shape[1] - 1) / 2)
    y, x = np.indices(shape, dtype=float)
    r = np.hypot(y - center[0], x - center[1])
    inside = 0.5 * (1.0 - np.tanh(2.0 * (r - radius) / width))
    return np.stack((1.0 - inside, inside))


def planar_interface(shape: tuple[int, int], width: float, angle: float = 0.0,
                     offset: float = 0.0) -> NDArray[np.float64]:
    y, x = np.indices(shape, dtype=float)
    yc, xc = (np.array(shape) - 1) / 2
    signed = (x - xc) * np.cos(angle) + (y - yc) * np.sin(angle) - offset
    phase = 0.5 * (1.0 + np.tanh(2.0 * signed / width))
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
    # Soft Voronoi indicators are normalized and exactly fill space.
    scale = max(width, 0.25) ** 2
    shifted = distance - distance.min(axis=0, keepdims=True)
    weights = np.exp(-shifted / scale)
    eta = weights / weights.sum(axis=0, keepdims=True)
    return eta, seeds, orientations


def equivalent_radius(phase: NDArray[np.float64], dx: float = 1.0) -> float:
    return float(np.sqrt(phase.sum() * dx**2 / np.pi))

