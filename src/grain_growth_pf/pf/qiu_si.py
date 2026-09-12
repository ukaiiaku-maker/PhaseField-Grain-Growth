"""Controlled anisotropic extension of the archived Qiu [100] SI update.

The A0 path reproduces the native nine-point capillary stencil and ordered-pair
algebra.  For anisotropic energy, the Cahn--Hoffman term is introduced as a
discrete variational correction which is identically zero for A0.  This keeps
native capillarity, elastic driving, and the antisymmetric ``eij`` term
separate and makes the scope of the extension auditable.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np
from numpy.typing import NDArray

from grain_growth_pf.pf.anisotropic import (
    _pair_homogeneous_norm_and_gradient,
    _pair_law,
)

Array = NDArray[np.float64]
A2_ANGULAR_SCALE = 1.0910512514090829


class QiuSIControl(str, Enum):
    ISO_NATIVE = "QIU_SI_4REF_ISO_NATIVE"
    A0_PORT = "QIU_SI_4REF_A0_PORT"
    ANISO_E = "QIU_SI_4REF_ANISO_E"
    ANISO_M = "QIU_SI_4REF_ANISO_M"
    ANISO_EM_INV = "QIU_SI_4REF_ANISO_EM_INV"


@dataclass(frozen=True)
class QiuSI4RefParameters:
    shape: tuple[int, int] = (500, 500)
    phases: int = 17
    dx: float = 1.0
    dy: float = 1.0
    dt: float = 0.1
    steps: int = 200_000
    interface_width: float = 5.0
    reference_count: int = 4
    output_cadence: int = 250
    eee: float = 20.0
    g_min: float = 0.65
    inclination_weight: float = 0.85
    support_power: int = 16
    mobility_exponent: float = 2.0

    @property
    def native_mobility(self) -> float:
        return 2.0 * np.pi**2 / (8.0 * self.interface_width)

    @property
    def reference_angle(self) -> float:
        return np.pi / self.reference_count


def qiu_native_beta(orientation_i: float, orientation_j: float) -> tuple[float, float]:
    """Exact scalar translation of the archived four-reference beta function."""
    theta = orientation_i - orientation_j
    if np.pi / 2 <= theta < np.pi:
        theta -= np.pi / 2
    if np.pi <= theta < 3 * np.pi / 2:
        theta -= np.pi
    if 3 * np.pi / 2 <= theta < 2 * np.pi:
        theta -= 3 * np.pi / 2
    if -np.pi < theta <= -np.pi / 2:
        theta += np.pi / 2
    if -3 * np.pi / 2 < theta <= -np.pi:
        theta += np.pi
    if -2 * np.pi < theta <= -3 * np.pi / 2:
        theta += 3 * np.pi / 2
    if abs(theta) <= np.deg2rad(36.7):
        first = 2.0 * np.tan(theta / 2.0)
    elif abs(theta) <= np.deg2rad(53.7):
        first = 0.0
    else:
        first = 0.0
    return float(first), float(-first)


def qiu_native_laplacian(field: Array, dx: float) -> Array:
    """Archived periodic nine-point stencil, including its single-dx scaling."""
    left = np.roll(field, 1, axis=0)
    right = np.roll(field, -1, axis=0)
    down = np.roll(field, 1, axis=1)
    up = np.roll(field, -1, axis=1)
    return (
        np.roll(left, 1, axis=1) + 4 * left + 4 * right
        + np.roll(right, -1, axis=1) + np.roll(left, -1, axis=1)
        + 4 * down + 4 * up + np.roll(right, 1, axis=1) - 20 * field
    ) / (6 * dx * dx)


def _local_presence(phi: Array) -> NDArray[np.bool_]:
    present = phi > 0.0
    return present | np.roll(present, 1, axis=1) | np.roll(present, -1, axis=1) | np.roll(present, 1, axis=2) | np.roll(present, -1, axis=2)


def qiu_native_capillary(phi_i: Array, phi_j: Array, lap_i: Array, lap_j: Array,
                         interface_width: float) -> Array:
    coefficient = np.pi**2 / (2.0 * interface_width**2)
    return phi_j * lap_i - phi_i * lap_j + coefficient * (phi_i - phi_j)


def qiu_pair_anisotropy_correction(
    phi_i: Array,
    phi_j: Array,
    orientation_i: float,
    orientation_j: float,
    dx: float,
    g_min: float,
    inclination_weight: float,
    support_power: int,
    angular_scale: float,
    energy_normalization: float,
) -> tuple[float, Array]:
    """Return ``E_A2-E_A0`` and the corresponding drive on phase ``i``.

    Forward bonds are counted once.  The returned drive is the negative
    derivative with respect to ``phi_i`` divided by cell area; its opposite is
    applied to ``phi_j``.  Support is held fixed during differentiation.
    """
    difference = phi_i - phi_j
    derivative = np.zeros_like(difference)
    energy = 0.0
    nx, ny = difference.shape
    for l in range(nx):
        lp = (l + 1) % nx
        for m in range(ny):
            mp = (m + 1) % ny
            px = (difference[lp, m] - difference[l, m]) / dx
            py = (difference[l, mp] - difference[l, m]) / dx
            norm, norm_x, norm_y = _pair_homogeneous_norm_and_gradient(
                px, py, orientation_i, orientation_j, 1.0, g_min,
                inclination_weight, support_power, angular_scale,
                energy_normalization,
            )
            qx = norm * norm_x
            qy = norm * norm_y
            energy += 0.25 * (norm * norm - px * px - py * py) * dx * dx
            flux_x = 0.5 * (qx - px) * dx
            flux_y = 0.5 * (qy - py) * dx
            derivative[l, m] -= flux_x + flux_y
            derivative[lp, m] += flux_x
            derivative[l, mp] += flux_y
    return float(energy), -derivative / (dx * dx)


def qiu_si_pairwise_rate(
    phi: Array,
    orientations: Array,
    signed_elastic_drive: Array,
    eij: Array,
    control: QiuSIControl,
    parameters: QiuSI4RefParameters = QiuSI4RefParameters(),
    *,
    angular_scale: float = A2_ANGULAR_SCALE,
    energy_normalization: float = 1.0,
    mobility_normalization: float = 1.0,
) -> tuple[Array, dict[str, float]]:
    """Compute the pre-renormalization Qiu increment as conservative pairs.

    ``signed_elastic_drive[i,j]`` is the already resolved native elastic PF
    force added to ordered pair ``(i,j)``.  It must remain antisymmetric; this
    function never recomputes, rescales, or duplicates it.
    """
    control = QiuSIControl(control)
    phi = np.asarray(phi, dtype=float)
    orientations = np.asarray(orientations, dtype=float)
    signed_elastic_drive = np.asarray(signed_elastic_drive, dtype=float)
    eij = np.asarray(eij, dtype=float)
    if np.any(~np.isfinite(phi)) or np.any(phi < 0.0):
        raise ValueError("phi must be finite and nonnegative before the Qiu update")
    phases, height, width = phi.shape
    if orientations.shape != (phases,) or eij.shape != (phases, phases):
        raise ValueError("orientation/eij dimensions do not match phi")
    if signed_elastic_drive.shape != (phases, phases, height, width):
        raise ValueError("signed_elastic_drive must have shape (phase, phase, y, x)")
    if not np.allclose(eij + eij.T, 0.0, rtol=0.0, atol=1e-14):
        raise ValueError("native eij must be antisymmetric")
    if not np.allclose(signed_elastic_drive + signed_elastic_drive.swapaxes(0, 1), 0.0, rtol=0.0, atol=1e-12):
        raise ValueError("native elastic pair drive must be antisymmetric")

    anisotropic_energy = control in (QiuSIControl.ANISO_E, QiuSIControl.ANISO_EM_INV)
    anisotropic_mobility = control in (QiuSIControl.ANISO_M, QiuSIControl.ANISO_EM_INV)
    if control is QiuSIControl.ISO_NATIVE:
        raise ValueError("ISO_NATIVE is executed only by the immutable archived source")
    laplacians = np.stack([qiu_native_laplacian(field, parameters.dx) for field in phi])
    local = _local_presence(phi)
    rate = np.zeros_like(phi)
    capillary_correction_energy = 0.0
    mobility_min = np.inf
    mobility_max = 0.0

    for i in range(phases - 1):
        for j in range(i + 1, phases):
            mask = local[i] & local[j]
            if not np.any(mask):
                continue
            capillary = qiu_native_capillary(
                phi[i], phi[j], laplacians[i], laplacians[j],
                parameters.interface_width,
            )
            if anisotropic_energy:
                correction_energy, correction = qiu_pair_anisotropy_correction(
                    phi[i], phi[j], orientations[i], orientations[j],
                    parameters.dx, parameters.g_min,
                    parameters.inclination_weight, parameters.support_power,
                    angular_scale, energy_normalization,
                )
                capillary = capillary + correction
                capillary_correction_energy += correction_energy

            mobility = np.full((height, width), parameters.native_mobility)
            if anisotropic_mobility:
                grad_ix = (np.roll(phi[i], -1, axis=0) - np.roll(phi[i], 1, axis=0)) / (2 * parameters.dx)
                grad_iy = (np.roll(phi[i], -1, axis=1) - np.roll(phi[i], 1, axis=1)) / (2 * parameters.dy)
                grad_jx = (np.roll(phi[j], -1, axis=0) - np.roll(phi[j], 1, axis=0)) / (2 * parameters.dx)
                grad_jy = (np.roll(phi[j], -1, axis=1) - np.roll(phi[j], 1, axis=1)) / (2 * parameters.dy)
                nx = phi[i] * grad_jx - phi[j] * grad_ix
                ny = phi[i] * grad_jy - phi[j] * grad_iy
                theta = np.arctan2(ny, nx)
                for y, x in zip(*np.nonzero(mask)):
                    mobility[y, x] = _pair_law(
                        theta[y, x], orientations[i], orientations[j], 1.0,
                        parameters.native_mobility, parameters.g_min,
                        parameters.inclination_weight, parameters.support_power,
                        parameters.mobility_exponent, angular_scale,
                        energy_normalization, mobility_normalization,
                    )[2]
            barrier = np.pi / parameters.interface_width * np.sqrt(phi[i] * phi[j]) * eij[i, j]
            complete_drive = capillary + signed_elastic_drive[i, j] + barrier
            exchange = np.where(mask, mobility * complete_drive, 0.0)
            rate[i] += exchange
            rate[j] -= exchange
            mobility_min = min(mobility_min, float(np.min(mobility[mask])))
            mobility_max = max(mobility_max, float(np.max(mobility[mask])))

    increment = parameters.dt * rate
    diagnostics = {
        "anisotropy_correction_energy": capillary_correction_energy,
        "phase_sum_rate_max_abs": float(np.max(np.abs(np.sum(rate, axis=0)))),
        "mobility_min": float(mobility_min if np.isfinite(mobility_min) else parameters.native_mobility),
        "mobility_max": float(mobility_max),
    }
    return increment, diagnostics
