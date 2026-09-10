from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def isotropic_lame_lambda(
    shear_modulus: float, poisson_ratio: float, constitutive_state: str
) -> float:
    """Return the in-plane Lamé coefficient for the selected 2-D closure."""
    if shear_modulus <= 0.0 or not (-1.0 < poisson_ratio < 0.5):
        raise ValueError("invalid isotropic elastic constants")
    if constitutive_state == "plane_strain":
        return 2.0 * shear_modulus * poisson_ratio / (1.0 - 2.0 * poisson_ratio)
    if constitutive_state == "plane_stress":
        return 2.0 * shear_modulus * poisson_ratio / (1.0 - poisson_ratio)
    raise ValueError("constitutive_state must be 'plane_stress' or 'plane_strain'")


class QiuFullFieldLegacy:
    """Historical periodic point-eigenstrain surrogate.

    This preserves the implementation used by the archived Phase-1 ``QIU``
    trajectory. It is not a faithful port of Qiu's current-geometry line
    construction and is retained under an explicit legacy name for replay.
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
        self.source_increment = np.zeros_like(self.eigenstrain)
        self.last_zero_mode_magnitude = 0.0
        self.last_equilibrium_residual = float("nan")

    def begin_source_step(self) -> None:
        """Reset the read-only accumulator for this coupling update."""
        self.source_increment.fill(0.0)

    def add_event(self, position: tuple[int, int], strain_increment: NDArray[np.float64]) -> None:
        inc = np.asarray(strain_increment, dtype=float)
        if inc.shape != (2, 2):
            raise ValueError("strain increment must be 2x2")
        y, x = position[0] % self.shape[0], position[1] % self.shape[1]
        symmetric = 0.5 * (inc + inc.T)
        self.eigenstrain[:, :, y, x] += symmetric
        self.source_increment[:, :, y, x] += symmetric

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
        self.last_zero_mode_magnitude = float(np.linalg.norm(sigma_hat[:, :, 0, 0]))
        sigma_hat[:, :, 0, 0] = 0.0
        divergence = np.empty((2, ny, nx), dtype=complex)
        for i in range(2):
            divergence[i] = yy * sigma_hat[i, 0] + xx * sigma_hat[i, 1]
        numerator = np.sqrt(np.sum(np.abs(divergence) ** 2, axis=0))
        sigma_norm = np.sqrt(np.sum(np.abs(sigma_hat) ** 2, axis=(0, 1)))
        denominator = np.sqrt(k2) * sigma_norm + np.finfo(float).eps
        nonzero = (xx != 0.0) | (yy != 0.0)
        self.last_equilibrium_residual = float(np.max(numerator[nonzero] / denominator[nonzero]))
        # Stress is the negative derivative of elastic energy with respect to
        # the imposed eigenstrain. The minus sign makes the self-field oppose
        # an additional like-signed transformation instead of amplifying it.
        self.stress = -np.fft.ifftn(sigma_hat, axes=(-2, -1)).real
        return self.stress

    def elastic_energy(self, grid_spacing: float = 1.0) -> float:
        """Historical surrogate energy implied by its self-stress sign."""
        return float(-0.5 * np.sum(self.stress * self.eigenstrain) * grid_spacing**2)

    def resolved_shear(self, position: tuple[int, int], tangent: NDArray[np.float64],
                       normal: NDArray[np.float64]) -> float:
        y, x = position[0] % self.shape[0], position[1] % self.shape[1]
        return float(np.asarray(tangent) @ self.stress[:, :, y, x] @ np.asarray(normal))


# Backward-compatible import for archived configs and reproducible replay.
QiuFullField = QiuFullFieldLegacy


class FFTEigenstrainV2:
    """Periodic isotropic microelasticity with exact Fourier equilibrium.

    Tensor index 0 is the array-row/y direction and index 1 is the
    array-column/x direction. At zero wavevector the default traction-free
    macroscopic closure sets the compatible mean strain equal to the imposed
    eigenstrain mean, hence the mean stress is zero.

    This is deliberately *not* called a Qiu reference implementation. Qiu's
    archived code reconstructs current-geometry line sources; this backend
    solves a distinct accumulated-eigenstrain model.
    """

    def __init__(
        self, shape: tuple[int, int], shear_modulus: float = 1.0,
        poisson_ratio: float = 0.3, *, grid_spacing: float = 1.0,
        constitutive_state: str = "plane_strain",
        zero_mode: str = "traction_free_mean_strain",
    ):
        if grid_spacing <= 0.0:
            raise ValueError("grid_spacing must be positive")
        if zero_mode not in {"traction_free_mean_strain", "fixed_mean_strain"}:
            raise ValueError("unknown zero_mode closure")
        self.shape = tuple(shape)
        self.mu = float(shear_modulus)
        self.nu = float(poisson_ratio)
        self.lam = isotropic_lame_lambda(self.mu, self.nu, constitutive_state)
        self.grid_spacing = float(grid_spacing)
        self.constitutive_state = constitutive_state
        self.zero_mode = zero_mode
        self.eigenstrain = np.zeros((2, 2, *self.shape), dtype=float)
        self.compatible_strain = np.zeros_like(self.eigenstrain)
        self.elastic_strain = np.zeros_like(self.eigenstrain)
        self.stress = np.zeros_like(self.eigenstrain)
        self.source_increment = np.zeros_like(self.eigenstrain)
        self.last_zero_mode_magnitude = 0.0
        self.last_equilibrium_residual = 0.0

    def begin_source_step(self) -> None:
        self.source_increment.fill(0.0)

    def add_event(
        self, position: tuple[int, int], strain_increment: NDArray[np.float64]
    ) -> None:
        increment = np.asarray(strain_increment, dtype=float)
        if increment.shape != (2, 2) or not np.all(np.isfinite(increment)):
            raise ValueError("strain increment must be a finite 2x2 tensor")
        y, x = position[0] % self.shape[0], position[1] % self.shape[1]
        symmetric = 0.5 * (increment + increment.T)
        self.eigenstrain[:, :, y, x] += symmetric
        self.source_increment[:, :, y, x] += symmetric

    def set_eigenstrain(self, eigenstrain: NDArray[np.float64]) -> None:
        value = np.asarray(eigenstrain, dtype=float)
        if value.shape != self.eigenstrain.shape or not np.all(np.isfinite(value)):
            raise ValueError("eigenstrain has the wrong shape or nonfinite values")
        self.eigenstrain = 0.5 * (value + np.swapaxes(value, 0, 1))

    def _hooke(self, strain_hat: NDArray[np.complex128]) -> NDArray[np.complex128]:
        stress_hat = 2.0 * self.mu * strain_hat.copy()
        trace = strain_hat[0, 0] + strain_hat[1, 1]
        stress_hat[0, 0] += self.lam * trace
        stress_hat[1, 1] += self.lam * trace
        return stress_hat

    def solve(self) -> NDArray[np.float64]:
        ny, nx = self.shape
        ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=self.grid_spacing)
        kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=self.grid_spacing)
        wave = np.stack(np.meshgrid(ky, kx, indexing="ij"))
        k2 = np.sum(wave * wave, axis=0)
        nonzero = k2 > 0.0
        eigen_hat = np.fft.fftn(self.eigenstrain, axes=(-2, -1))

        # b_i = k_j C_ijkl epsilon*_kl.
        trace_star = eigen_hat[0, 0] + eigen_hat[1, 1]
        b = np.empty((2, ny, nx), dtype=complex)
        for i in range(2):
            b[i] = self.lam * wave[i] * trace_star
            for j in range(2):
                b[i] += 2.0 * self.mu * wave[j] * eigen_hat[i, j]

        # A^-1 = (mu k^2)^-1 P_T + ((lambda+2mu) k^2)^-1 P_L.
        inverse_b = np.zeros_like(b)
        longitudinal_dot = np.sum(wave * b, axis=0)
        for i in range(2):
            transverse = b[i].copy()
            transverse[nonzero] -= (
                wave[i, nonzero] * longitudinal_dot[nonzero] / k2[nonzero]
            )
            inverse_b[i, nonzero] = (
                transverse[nonzero] / (self.mu * k2[nonzero])
                + wave[i, nonzero] * longitudinal_dot[nonzero]
                / ((self.lam + 2.0 * self.mu) * k2[nonzero] ** 2)
            )

        compatible_hat = np.zeros_like(eigen_hat)
        for i in range(2):
            for j in range(2):
                compatible_hat[i, j] = 0.5 * (
                    wave[i] * inverse_b[j] + wave[j] * inverse_b[i]
                )
        unconstrained_zero_stress = self._hooke(-eigen_hat[:, :, 0:1, 0:1])
        self.last_zero_mode_magnitude = float(np.linalg.norm(unconstrained_zero_stress))
        if self.zero_mode == "traction_free_mean_strain":
            compatible_hat[:, :, 0, 0] = eigen_hat[:, :, 0, 0]

        elastic_hat = compatible_hat - eigen_hat
        stress_hat = self._hooke(elastic_hat)
        divergence = np.einsum("jYX,ijYX->iYX", wave, stress_hat)
        numerator = np.sqrt(np.sum(np.abs(divergence) ** 2, axis=0))
        stress_norm = np.sqrt(np.sum(np.abs(stress_hat) ** 2, axis=(0, 1)))
        denominator = np.sqrt(k2) * stress_norm + np.finfo(float).eps
        self.last_equilibrium_residual = float(
            np.max(numerator[nonzero] / denominator[nonzero])
        )
        self.compatible_strain = np.fft.ifftn(
            compatible_hat, axes=(-2, -1)
        ).real
        self.elastic_strain = np.fft.ifftn(elastic_hat, axes=(-2, -1)).real
        self.stress = np.fft.ifftn(stress_hat, axes=(-2, -1)).real
        return self.stress

    def elastic_energy(self, grid_spacing: float | None = None) -> float:
        dx = self.grid_spacing if grid_spacing is None else float(grid_spacing)
        return float(0.5 * np.sum(self.elastic_strain * self.stress) * dx**2)

    def resolved_shear(
        self, position: tuple[int, int], tangent: NDArray[np.float64],
        normal: NDArray[np.float64],
    ) -> float:
        y, x = position[0] % self.shape[0], position[1] % self.shape[1]
        return float(np.asarray(tangent) @ self.stress[:, :, y, x] @ np.asarray(normal))
