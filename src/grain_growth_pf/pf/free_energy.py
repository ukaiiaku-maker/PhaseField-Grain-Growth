from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]


def laplacian(field: Array, dx: float, boundary: str = "periodic") -> Array:
    """Second-order five-point Laplacian."""
    if boundary == "periodic":
        return (
            np.roll(field, 1, -2) + np.roll(field, -1, -2)
            + np.roll(field, 1, -1) + np.roll(field, -1, -1) - 4.0 * field
        ) / dx**2
    padded = np.pad(field, ((0, 0), (1, 1), (1, 1)), mode="edge")
    return (
        padded[:, 2:, 1:-1] + padded[:, :-2, 1:-1]
        + padded[:, 1:-1, 2:] + padded[:, 1:-1, :-2]
        - 4.0 * field
    ) / dx**2


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
                stored_energy: float = 0.0) -> float:
    kappa, well = coefficients(gamma, width)
    gx = (np.roll(eta, -1, axis=1) - eta) / dx
    gy = (np.roll(eta, -1, axis=2) - eta) / dx
    density = 0.5 * kappa * (gx * gx + gy * gy) + well * eta**2 * (1.0 - eta)**2
    # The sum over phases counts both sides of an interface. The 1/2 matches
    # the pair-interface convention used by the dynamics.
    return float(0.5 * density.sum() * dx**2 + stored_energy)

