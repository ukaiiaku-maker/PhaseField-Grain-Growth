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
    partner_phase: NDArray[np.float64] | None = None,
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

    # A radius-two stencil suppresses the one-cell lattice oscillation that
    # otherwise dominates a direct divergence-of-normal estimate.
    radius = 2
    yc = index(y, height)
    xc = index(x, width)
    yp = index(y + radius, height)
    ym = index(y - radius, height)
    xp = index(x + radius, width)
    xm = index(x - radius, width)
    span = 2.0 * radius * dx
    gradient_y = (field[yp, xc] - field[ym, xc]) / span
    gradient_x = (field[yc, xp] - field[yc, xm]) / span
    second_scale = (radius * dx) ** 2
    second_y = (field[yp, xc] - 2.0 * field[yc, xc] + field[ym, xc]) / second_scale
    second_x = (field[yc, xp] - 2.0 * field[yc, xc] + field[yc, xm]) / second_scale
    mixed = (
        field[yp, xp] - field[yp, xm] - field[ym, xp] + field[ym, xm]
    ) / (4.0 * second_scale)
    gradient_magnitude = np.hypot(gradient_y, gradient_x)
    curvature = (
        second_y * gradient_x**2
        - 2.0 * mixed * gradient_y * gradient_x
        + second_x * gradient_y**2
    ) / np.maximum(gradient_magnitude**3, 1e-14)
    normal_y = -gradient_y / np.maximum(gradient_magnitude, 1e-14)
    normal_x = -gradient_x / np.maximum(gradient_magnitude, 1e-14)
    valid = gradient_magnitude > 1e-10
    valid &= (field[yc, xc] > 0.05) & (field[yc, xc] < 0.95)
    if partner_phase is not None:
        partner = np.asarray(partner_phase, float)
        if partner.shape != field.shape:
            raise ValueError("partner phase has the wrong shape")
        valid &= field[yc, xc] + partner[yc, xc] > 0.90
    if not np.any(valid):
        return 0.0, 0.0, (0.0, 0.0)
    yp_one = index(y + 1, height)
    ym_one = index(y - 1, height)
    xp_one = index(x + 1, width)
    xm_one = index(x - 1, width)
    velocity_gradient = np.hypot(
        (field[yp_one, xc] - field[ym_one, xc]) / (2.0 * dx),
        (field[yc, xp_one] - field[yc, xm_one]) / (2.0 * dx),
    )
    velocity = (field[yc, xc] - previous[yc, xc]) / elapsed
    velocity /= np.maximum(velocity_gradient, 1e-14)
    mean_normal = np.array([np.mean(normal_y[valid]), np.mean(normal_x[valid])])
    normal_norm = np.linalg.norm(mean_normal)
    if normal_norm > 1e-6:
        mean_normal /= normal_norm
    else:
        mean_normal[:] = 0.0
    return (
        float(np.mean(curvature[valid])),
        float(np.mean(velocity[valid])),
        (float(mean_normal[0]), float(mean_normal[1])),
    )
