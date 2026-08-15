from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def interface_kinematics(
    phase: NDArray[np.float64],
    previous_phase: NDArray[np.float64],
    points: NDArray[np.float64],
    elapsed: float,
    dx: float,
    periodic: bool = True,
) -> tuple[float, float, tuple[float, float]]:
    """Measure signed local curvature and velocity from a diffuse level set.

    The normal points out of ``phase``. Curvature is signed so intrinsic
    capillary motion follows ``v = M gamma kappa``: a convex shrinking grain
    has both negative curvature and negative outward velocity.
    """
    field = np.asarray(phase, float)
    previous = np.asarray(previous_phase, float)
    coordinates = np.asarray(points, float).astype(int)
    if field.ndim != 2 or previous.shape != field.shape or coordinates.ndim != 2:
        raise ValueError("invalid phase-field kinematics inputs")
    if elapsed <= 0 or dx <= 0 or len(coordinates) == 0:
        return 0.0, 0.0, (0.0, 0.0)
    height, width = field.shape
    y = coordinates[:, 0]
    x = coordinates[:, 1]

    def index(values: NDArray[np.int64], size: int) -> NDArray[np.int64]:
        return values % size if periodic else np.clip(values, 0, size - 1)

    def outward_normal(y_offset: int, x_offset: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        yc = index(y + y_offset, height)
        xc = index(x + x_offset, width)
        yp = index(y + y_offset + 1, height)
        ym = index(y + y_offset - 1, height)
        xp = index(x + x_offset + 1, width)
        xm = index(x + x_offset - 1, width)
        gradient_y = (field[yp, xc] - field[ym, xc]) / (2.0 * dx)
        gradient_x = (field[yc, xp] - field[yc, xm]) / (2.0 * dx)
        magnitude = np.hypot(gradient_y, gradient_x)
        safe = np.maximum(magnitude, 1e-14)
        return -gradient_y / safe, -gradient_x / safe, magnitude

    normal_y, normal_x, gradient_magnitude = outward_normal(0, 0)
    normal_y_plus, _, _ = outward_normal(1, 0)
    normal_y_minus, _, _ = outward_normal(-1, 0)
    _, normal_x_plus, _ = outward_normal(0, 1)
    _, normal_x_minus, _ = outward_normal(0, -1)
    divergence = (
        normal_y_plus - normal_y_minus + normal_x_plus - normal_x_minus
    ) / (2.0 * dx)
    valid = gradient_magnitude > 1e-10
    if not np.any(valid):
        return 0.0, 0.0, (0.0, 0.0)
    yc = index(y, height)
    xc = index(x, width)
    velocity = (field[yc, xc] - previous[yc, xc]) / elapsed
    velocity /= np.maximum(gradient_magnitude, 1e-14)
    mean_normal = np.array([np.mean(normal_y[valid]), np.mean(normal_x[valid])])
    normal_norm = np.linalg.norm(mean_normal)
    if normal_norm > 1e-6:
        mean_normal /= normal_norm
    else:
        mean_normal[:] = 0.0
    return (
        float(np.mean(-divergence[valid])),
        float(np.mean(velocity[valid])),
        (float(mean_normal[0]), float(mean_normal[1])),
    )
