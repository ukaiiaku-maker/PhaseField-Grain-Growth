from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]


def laplacian(field: Array, dx: float, boundary: str = "periodic") -> Array:
    """Nine-point, second-order isotropic Laplacian used by Qiu et al."""
    if boundary == "periodic":
        cardinal = (
            np.roll(field, 1, -2) + np.roll(field, -1, -2)
            + np.roll(field, 1, -1) + np.roll(field, -1, -1)
        )
        diagonal = (
            np.roll(np.roll(field, 1, -2), 1, -1)
            + np.roll(np.roll(field, 1, -2), -1, -1)
            + np.roll(np.roll(field, -1, -2), 1, -1)
            + np.roll(np.roll(field, -1, -2), -1, -1)
        )
        return (4.0 * cardinal + diagonal - 20.0 * field) / (6.0 * dx**2)
    padded = np.pad(field, ((0, 0), (1, 1), (1, 1)), mode="edge")
    cardinal = (
        padded[:, 2:, 1:-1] + padded[:, :-2, 1:-1]
        + padded[:, 1:-1, 2:] + padded[:, 1:-1, :-2]
    )
    diagonal = (
        padded[:, 2:, 2:] + padded[:, 2:, :-2]
        + padded[:, :-2, 2:] + padded[:, :-2, :-2]
    )
    return (4.0 * cardinal + diagonal - 20.0 * field) / (6.0 * dx**2)


def coefficients(gamma: float, width: float) -> tuple[float, float]:
    """Return gradient and double-well coefficients with surface energy gamma.

    For f = kappa |grad eta|^2/2 + W eta^2(1-eta)^2, the isolated
    interface energy is sqrt(2*kappa*W)/6. The chosen ratio gives a diffuse
    width proportional to ``width`` while preserving that energy exactly.
    """
    kappa = 3.0 * gamma * width
    well = 6.0 * gamma / width
    return kappa, well


def chemical_potential(eta: Array, gamma: float, width: float, dx: float,
                       boundary: str = "periodic") -> Array:
    kappa, well = coefficients(gamma, width)
    bulk = 2.0 * well * eta * (1.0 - eta) * (1.0 - 2.0 * eta)
    return bulk - kappa * laplacian(eta, dx, boundary)


def free_energy(eta: Array, gamma: float, width: float, dx: float,
                stored_energy: float = 0.0,
                boundary: str = "periodic") -> float:
    """Pairwise double-obstacle interfacial energy.

    The efficient pair sums are algebraically equivalent to summing over
    ``i < j``.  The normalization gives an isolated equilibrium interface
    energy ``gamma`` for the compact sinusoidal profile of total width
    ``width``.
    """
    if boundary == "periodic":
        gx = (np.roll(eta, -1, axis=1) - eta) / dx
        gy = (np.roll(eta, -1, axis=2) - eta) / dx
    else:
        gx = np.zeros_like(eta)
        gy = np.zeros_like(eta)
        gx[:, :-1] = (eta[:, 1:] - eta[:, :-1]) / dx
        gy[:, :, :-1] = (eta[:, :, 1:] - eta[:, :, :-1]) / dx
    phase_sum = eta.sum(axis=0)
    pair_potential = 0.5 * (phase_sum**2 - np.sum(eta**2, axis=0))
    pair_gradient = 0.5 * (
        gx.sum(axis=0) ** 2 + gy.sum(axis=0) ** 2
        - np.sum(gx**2 + gy**2, axis=0)
    )
    density = (4.0 * gamma / width) * (
        pair_potential - width**2 / np.pi**2 * pair_gradient
    )
    return float(density.sum() * dx**2 + stored_energy)
