from __future__ import annotations

import numpy as np


def event_shear_increment(burgers_parallel: float, swept_length: float, rve_area: float) -> float:
    if rve_area <= 0:
        raise ValueError("RVE area must be positive")
    return float(burgers_parallel * swept_length / rve_area)


def infer_rotation(tangential_displacements: np.ndarray, positions: np.ndarray,
                   centroid: tuple[float, float]) -> float:
    """Least-squares small rotation from spatial event displacements."""
    r = np.asarray(positions, dtype=float) - np.asarray(centroid, dtype=float)
    u = np.asarray(tangential_displacements, dtype=float)
    denom = float(np.sum(r * r))
    return 0.0 if denom == 0 else float(np.sum(r[:, 0] * u[:, 1] - r[:, 1] * u[:, 0]) / denom)

