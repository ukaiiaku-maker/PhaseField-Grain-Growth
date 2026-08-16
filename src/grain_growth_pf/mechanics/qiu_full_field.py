from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


class QiuFullField:
    """Periodic nonlocal elastic backend from event eigenstrain.

    This is an independent Fourier-space isotropic Green-operator closure. It
    serves the same scientific role as Qiu's real-space line-disconnection
    superposition while accepting event-resolved plastic eigenstrain fields.
    The k=0 mode is removed (traction-free mean strain).
    """

    def __init__(self, shape: tuple[int, int], shear_modulus: float = 1.0,
                 poisson_ratio: float = 0.3):
        if not (-1.0 < poisson_ratio < 0.5) or shear_modulus <= 0:
            raise ValueError("invalid isotropic elastic constants")
        self.shape = shape
        self.mu = shear_modulus
        self.nu = poisson_ratio
        self.eigenstrain = np.zeros((2, 2, *shape), dtype=float)
        self.stress = np.zeros_like(self.eigenstrain)

    def add_event(self, position: tuple[int, int], strain_increment: NDArray[np.float64]) -> None:
        inc = np.asarray(strain_increment, dtype=float)
        if inc.shape != (2, 2):
            raise ValueError("strain increment must be 2x2")
        y, x = position[0] % self.shape[0], position[1] % self.shape[1]
        self.eigenstrain[:, :, y, x] += 0.5 * (inc + inc.T)

    def solve(self) -> NDArray[np.float64]:
        ny, nx = self.shape
        ky = 2 * np.pi * np.fft.fftfreq(ny)
        kx = 2 * np.pi * np.fft.fftfreq(nx)
        yy, xx = np.meshgrid(ky, kx, indexing="ij")
        k2 = xx * xx + yy * yy
        k2[0, 0] = 1.0
        # Incompatibility projection: retain the nonlocal, divergence-balanced
        # part of eigenstrain and apply plane-strain Hooke response.
        eps = np.fft.fftn(self.eigenstrain, axes=(-2, -1))
        n = np.stack((yy / np.sqrt(k2), xx / np.sqrt(k2)))
        projected = np.empty_like(eps)
        identity = np.eye(2)
        for i in range(2):
            for j in range(2):
                projected[i, j] = eps[i, j]
                for a in range(2):
                    projected[i, j] -= 0.5 * (
                        n[i] * n[a] * eps[a, j] + n[j] * n[a] * eps[i, a]
                    )
        trace = projected[0, 0] + projected[1, 1]
        lam = 2 * self.mu * self.nu / (1 - 2 * self.nu)
        sigma_hat = 2 * self.mu * projected
        sigma_hat[0, 0] += lam * trace
        sigma_hat[1, 1] += lam * trace
        sigma_hat[:, :, 0, 0] = 0.0
        # Stress is the negative derivative of elastic energy with respect to
        # the imposed eigenstrain. The minus sign makes the self-field oppose
        # an additional like-signed transformation instead of amplifying it.
        self.stress = -np.fft.ifftn(sigma_hat, axes=(-2, -1)).real
        return self.stress

    def resolved_shear(self, position: tuple[int, int], tangent: NDArray[np.float64],
                       normal: NDArray[np.float64]) -> float:
        y, x = position[0] % self.shape[0], position[1] % self.shape[1]
        return float(np.asarray(tangent) @ self.stress[:, :, y, x] @ np.asarray(normal))
