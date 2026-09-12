from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from numpy.typing import NDArray

from grain_growth_pf.config import PFConfig
from grain_growth_pf.mechanics.anisotropy import (
    LADDER,
    angular_normalization,
    support_derivatives,
)
from .anisotropic import anisotropic_energy_gradient, anisotropic_pairwise_step
from .free_energy import free_energy
from .kernels import pairwise_obstacle_step

Array = NDArray[np.float64]
DrivingCallback = Callable[[Array, float], Array]


def project_simplex(values: Array) -> Array:
    """Project each phase vector onto the probability simplex exactly."""
    moved = np.moveaxis(values, 0, -1)
    flat = moved.reshape(-1, moved.shape[-1])
    u = np.sort(flat, axis=1)[:, ::-1]
    cssv = np.cumsum(u, axis=1) - 1.0
    ind = np.arange(1, u.shape[1] + 1)
    cond = u - cssv / ind > 0
    rho = cond.sum(axis=1) - 1
    theta = cssv[np.arange(len(flat)), rho] / (rho + 1)
    projected = np.maximum(flat - theta[:, None], 0.0)
    return np.moveaxis(projected.reshape(moved.shape), -1, 0)


@dataclass
class StepDiagnostics:
    time: float
    step: int
    dt: float
    interfacial_energy: float
    max_constraint_error: float


class MultiphaseFieldSolver:
    """Pairwise constrained multiphase-field solver following Qiu et al.

    The capillary force is the published pairwise double-obstacle form.  Only
    locally present phases and their one-stencil-cell halos participate, which
    permits neighbor switching without remote phase nucleation.
    """

    def __init__(self, eta: Array, config: PFConfig,
                 driving: DrivingCallback | None = None,
                 orientations: Array | None = None):
        eta = np.asarray(eta, dtype=float)
        if eta.ndim != 3 or eta.shape[1:] != config.shape:
            raise ValueError("eta must have shape (n_grains, *config.shape)")
        sums = eta.sum(axis=0)
        already_on_simplex = bool(
            np.all(np.isfinite(eta))
            and np.min(eta) >= 0.0
            and np.max(np.abs(sums - 1.0)) <= 1e-12
        )
        self.eta = eta.copy() if already_on_simplex else project_simplex(eta)
        self.config = config
        self.driving = driving
        self.orientations = (
            None if orientations is None else np.asarray(orientations, dtype=float).copy()
        )
        self.anisotropic = bool(
            config.anisotropy_strength not in {None, "A0_ISOTROPIC"}
            and (config.anisotropic_energy or config.anisotropic_mobility)
        )
        phase_maxima = np.max(self.eta, axis=(1, 2))
        self.active_phases = (
            phase_maxima > 0.0 if self.anisotropic
            else phase_maxima >= config.grain_extinction_threshold
        )
        if not np.any(self.active_phases):
            raise ValueError("initial condition contains no active grain")
        if self.anisotropic and (
            self.orientations is None or self.orientations.shape != (len(self.eta),)
        ):
            raise ValueError("anisotropic PF evolution requires one orientation per phase")
        self._last_capillary_potential: Array | None = None
        self._last_pre_step_energy = float("nan")
        self.mobility_scale = np.ones(config.shape, dtype=float)
        self.time = 0.0
        self.step_number = 0

    @property
    def labels(self) -> NDArray[np.int64]:
        return np.argmax(self.eta, axis=0)

    def stable_dt(self) -> float:
        # Explicit diffusion stability bound in 2-D; the factor 0.18 leaves
        # margin for the local double-well term.
        kappa = 3.0 * self.config.gb_energy * self.config.interface_width
        kinetic = self.config.intrinsic_mobility / (3.0 * self.config.interface_width)
        base = 0.18 * self.config.grid_spacing**2 / max(
            kinetic * kappa, np.finfo(float).tiny
        )
        if not self.anisotropic:
            return base
        strength = LADDER[str(self.config.anisotropy_strength)]
        angles = np.arange(4096, dtype=float) * (2.0 * np.pi / 4096.0)
        support_stiffness = support_derivatives(
            angles, strength.support_power
        )[3]
        angular_scale = angular_normalization(strength)
        gamma_min = (
            self.config.gb_energy
            * self.config.anisotropy_energy_normalization
            * strength.g_min
            * angular_scale
            * (1.0 - strength.inclination_weight)
        )
        stiffness_max = (
            self.config.gb_energy
            * self.config.anisotropy_energy_normalization
            * angular_scale
            * (
                1.0 - strength.inclination_weight
                + strength.inclination_weight * float(np.max(support_stiffness))
            )
        )
        mobility_max = (
            self.config.intrinsic_mobility
            * self.config.anisotropy_mobility_normalization
            * (gamma_min / self.config.gb_energy) ** (-strength.mobility_exponent)
            if self.config.anisotropic_mobility else self.config.intrinsic_mobility
        )
        capillary_max = stiffness_max if self.config.anisotropic_energy else self.config.gb_energy
        return base / max(
            1.0,
            mobility_max * capillary_max
            / (self.config.intrinsic_mobility * self.config.gb_energy),
        )

    def step(self, dt: float | None = None, *, compute_energy: bool = True) -> StepDiagnostics:
        cfg = self.config
        requested = cfg.time_step if dt is None else dt
        used_dt = min(requested, self.stable_dt()) if cfg.adaptive_stepping else requested
        external = np.empty((1, 1, 1), dtype=float)
        use_external = self.driving is not None
        if self.driving is not None:
            external = np.asarray(self.driving(self.eta, self.time), dtype=float)
            if external.shape != self.eta.shape:
                raise ValueError("driving callback returned the wrong shape")
        if self.anisotropic:
            strength = LADDER[str(cfg.anisotropy_strength)]
            self.eta, self._last_pre_step_energy, self._last_capillary_potential = anisotropic_pairwise_step(
                self.eta, self.active_phases, self.orientations,
                self.mobility_scale, external, use_external, used_dt,
                cfg.gb_energy, cfg.intrinsic_mobility, cfg.interface_width,
                cfg.grid_spacing, cfg.boundary_conditions == "periodic",
                strength.g_min, strength.inclination_weight,
                strength.support_power, strength.mobility_exponent,
                angular_normalization(strength),
                cfg.anisotropy_energy_normalization,
                cfg.anisotropy_mobility_normalization,
                cfg.anisotropic_energy, cfg.anisotropic_mobility,
            )
        else:
            self.eta = pairwise_obstacle_step(
                self.eta,
                self.active_phases,
                self.mobility_scale,
                external,
                use_external,
                used_dt,
                cfg.intrinsic_mobility,
                cfg.gb_energy,
                cfg.interface_width,
                cfg.grid_spacing,
                cfg.boundary_conditions == "periodic",
            )
        if np.any(~np.isfinite(self.eta)):
            raise FloatingPointError("PF update produced a nonfinite phase value")
        phase_maxima = np.max(self.eta, axis=(1, 2))
        extinct = self.active_phases & (
            phase_maxima == 0.0 if self.anisotropic
            else phase_maxima < cfg.grain_extinction_threshold
        )
        if np.any(extinct) and np.count_nonzero(self.active_phases) > np.count_nonzero(extinct):
            self.active_phases[extinct] = False
            if not self.anisotropic:
                self.eta[extinct] = 0.0
                self.eta /= self.eta.sum(axis=0, keepdims=True)
        self.time += used_dt
        self.step_number += 1
        return StepDiagnostics(
            self.time, self.step_number, used_dt,
            (
                self._anisotropic_energy()
                if compute_energy and self.anisotropic
                else free_energy(
                    self.eta, cfg.gb_energy, cfg.interface_width, cfg.grid_spacing,
                    boundary=cfg.boundary_conditions,
                )
                if compute_energy else float("nan")
            ),
            float(np.max(np.abs(self.eta.sum(axis=0) - 1.0))),
        )

    def _anisotropic_energy(self) -> float:
        cfg = self.config
        strength = LADDER[str(cfg.anisotropy_strength)]
        energy, _ = anisotropic_energy_gradient(
            self.eta, self.active_phases, self.orientations,
            cfg.gb_energy, cfg.intrinsic_mobility, cfg.interface_width,
            cfg.grid_spacing, cfg.boundary_conditions == "periodic",
            strength.g_min, strength.inclination_weight,
            strength.support_power, strength.mobility_exponent,
            angular_normalization(strength),
            cfg.anisotropy_energy_normalization,
            cfg.anisotropy_mobility_normalization,
            cfg.anisotropic_energy,
        )
        return float(energy)

    def capillary_pressure(
        self, grain_i: int, grain_j: int, points: NDArray[np.integer]
    ) -> float | None:
        """Return the pair thermodynamic drive used by the last PF update."""
        if not self.anisotropic or self._last_capillary_potential is None:
            return None
        sample = np.asarray(points, dtype=int)
        if sample.ndim != 2 or sample.shape[1] != 2 or not len(sample):
            return None
        y = sample[:, 0] % self.config.shape[0]
        x = sample[:, 1] % self.config.shape[1]
        scale = np.pi**2 / (
            4.0 * self.config.interface_width * self.config.grid_spacing**2
        )
        values = scale * (
            self._last_capillary_potential[grain_j, y, x]
            - self._last_capillary_potential[grain_i, y, x]
        )
        finite = values[np.isfinite(values)]
        return float(np.mean(finite)) if finite.size else None

    def run(self, steps: int, callback: Callable[["MultiphaseFieldSolver", StepDiagnostics], None] | None = None) -> list[StepDiagnostics]:
        records: list[StepDiagnostics] = []
        for _ in range(steps):
            diag = self.step()
            records.append(diag)
            if callback is not None:
                callback(self, diag)
        return records

    def set_mobility_scale(self, scale: Array | float) -> None:
        value = np.asarray(scale, dtype=float)
        if value.ndim == 0:
            value = np.full(self.config.shape, float(value))
        if value.shape != self.config.shape or np.any(value < 0) or np.any(~np.isfinite(value)):
            raise ValueError("mobility scale must be a finite nonnegative spatial field")
        self.mobility_scale = value.copy()

    def state_dict(self) -> dict[str, object]:
        return {"eta": self.eta.copy(), "time": self.time, "step_number": self.step_number,
                "mobility_scale": self.mobility_scale.copy(),
                "active_phases": self.active_phases.copy(),
                "orientations": (
                    None if self.orientations is None else self.orientations.copy()
                )}

    def load_state_dict(self, state: dict[str, object]) -> None:
        self.eta = np.asarray(state["eta"], dtype=float).copy()
        self.time = float(state["time"])
        self.step_number = int(state["step_number"])
        self.mobility_scale = np.asarray(state.get("mobility_scale", np.ones(self.config.shape)), dtype=float).copy()
        phase_maxima = np.max(self.eta, axis=(1, 2))
        default_active = (
            phase_maxima > 0.0 if self.anisotropic
            else phase_maxima >= self.config.grain_extinction_threshold
        )
        self.active_phases = np.asarray(
            state.get("active_phases", default_active), dtype=bool
        ).copy()
        restored_orientations = state.get("orientations")
        if restored_orientations is not None:
            self.orientations = np.asarray(restored_orientations, dtype=float).copy()
