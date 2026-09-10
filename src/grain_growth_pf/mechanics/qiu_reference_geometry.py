from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def two_reference_coupling_factors(
    orientation_1: float, orientation_2: float
) -> tuple[float, float]:
    """Exact scalar ``beta`` law in archived Polycrystals/functions_2ref.py.

    Orientations are radians. The directed difference and endpoint handling are
    intentionally retained because the two returned reference coefficients are
    ordered, not an unsigned grain-boundary property.
    """
    theta = float(orientation_1 - orientation_2)
    if theta >= np.pi:
        theta -= np.pi
    elif theta <= -np.pi:
        theta += np.pi

    def branch(value: float) -> float:
        magnitude = abs(value)
        if magnitude <= np.deg2rad(70.0):
            return float(2.0 * np.tan(value / 2.0))
        if magnitude <= np.deg2rad(150.0):
            return float(2.0 * np.tan(
                value / 2.0 - np.sign(value) * np.deg2rad(109.38) / 2.0
            ))
        return float(-2.0 * np.tan(np.sign(value) * np.pi / 2.0 - value / 2.0))

    complementary = np.sign(theta) * np.pi - theta
    return 0.0333 * branch(theta), 0.0333 * branch(complementary)


@dataclass(frozen=True)
class TwoReferenceDecomposition:
    reference_angle_1: float
    reference_angle_2: float
    coupling_1: float
    coupling_2: float


def two_reference_decomposition(
    orientation_1: float, orientation_2: float,
    pair_normal_yx: np.ndarray,
) -> TwoReferenceDecomposition:
    """Reproduce Qiu's two-reference choice for a local pair normal.

    Inputs and returned reference angles use the repository's y/x tensor order.
    The archived source calls its first array coordinate x; this conversion is
    kept at this API boundary and tested explicitly.
    """
    normal = np.asarray(pair_normal_yx, dtype=float)
    if normal.shape != (2,) or not np.all(np.isfinite(normal)):
        raise ValueError("pair_normal_yx must be a finite two-vector")
    mean_orientation = 0.5 * (orientation_1 + orientation_2)
    # Archived dnormal_x/dnormal_y after mapping repository (y,x) to its (x,y).
    archived_x, archived_y = normal[1], normal[0]
    rotated_x = np.cos(mean_orientation) * archived_x + np.sin(mean_orientation) * archived_y
    rotated_y = -np.sin(mean_orientation) * archived_x + np.cos(mean_orientation) * archived_y
    beta_1, beta_2 = two_reference_coupling_factors(orientation_1, orientation_2)
    first_sector = rotated_y != 0.0 and -rotated_x / rotated_y >= 0.0
    if first_sector:
        return TwoReferenceDecomposition(
            mean_orientation, mean_orientation + np.pi / 2.0, beta_1, beta_2
        )
    return TwoReferenceDecomposition(
        mean_orientation + np.pi / 2.0, mean_orientation, beta_2, beta_1
    )


def resolved_reference_shear(stress_yx: np.ndarray, angle_from_archived_x: float) -> float:
    """Evaluate the archived resolved-shear expression with explicit axes."""
    stress = np.asarray(stress_yx, dtype=float)
    if stress.shape != (2, 2):
        raise ValueError("stress_yx must be a 2x2 tensor")
    # archived sigma11=xx, sigma22=yy, sigma12=xy
    sigma_xx = stress[1, 1]
    sigma_yy = stress[0, 0]
    sigma_xy = 0.5 * (stress[0, 1] + stress[1, 0])
    return float(
        0.5 * (sigma_yy - sigma_xx) * np.sin(2.0 * angle_from_archived_x)
        + sigma_xy * np.cos(2.0 * angle_from_archived_x)
    )


def qiu_reference_elastic_projection(
    stress_yx: np.ndarray, decomposition: TwoReferenceDecomposition
) -> float:
    """The stress/coupling factor in archived ``E_elastic`` (sans PF weights)."""
    return float(-(
        resolved_reference_shear(stress_yx, decomposition.reference_angle_1)
        * decomposition.coupling_1
        + resolved_reference_shear(stress_yx, decomposition.reference_angle_2)
        * decomposition.coupling_2
    ))
