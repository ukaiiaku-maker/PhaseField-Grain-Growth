from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from numpy.typing import NDArray

from grain_growth_pf.config import PFConfig
from .free_energy import free_energy, laplacian

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
                 driving: DrivingCallback | None = None):
        eta = np.asarray(eta, dtype=float)
        if eta.ndim != 3 or eta.shape[1:] != config.shape:
            raise ValueError("eta must have shape (n_grains, *config.shape)")
        self.eta = project_simplex(eta)
        self.active_phases = np.max(self.eta, axis=(1, 2)) >= config.grain_extinction_threshold
        if not np.any(self.active_phases):
            raise ValueError("initial condition contains no active grain")
        self.config = config
        self.driving = driving
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
        return 0.18 * self.config.grid_spacing**2 / max(
            kinetic * kappa, np.finfo(float).tiny
        )

    def step(self, dt: float | None = None) -> StepDiagnostics:
        cfg = self.config
        requested = cfg.time_step if dt is None else dt
        used_dt = min(requested, self.stable_dt()) if cfg.adaptive_stepping else requested
        active_indices = np.flatnonzero(self.active_phases)
        eta_active = self.eta[active_indices]
        lap = laplacian(eta_active, cfg.grid_spacing, cfg.boundary_conditions)
        present = eta_active > 1e-14
        if cfg.boundary_conditions == "periodic":
            support = (
                present | np.roll(present, 1, -2) | np.roll(present, -1, -2)
                | np.roll(present, 1, -1) | np.roll(present, -1, -1)
            )
        else:
            padded = np.pad(present, ((0, 0), (1, 1), (1, 1)), mode="constant")
            support = (
                present | padded[:, 2:, 1:-1] | padded[:, :-2, 1:-1]
                | padded[:, 1:-1, 2:] | padded[:, 1:-1, :-2]
            )
        count = np.maximum(support.sum(axis=0, keepdims=True), 1)
        phase_sum = (eta_active * support).sum(axis=0, keepdims=True)
        lap_sum = (lap * support).sum(axis=0, keepdims=True)
        # Match the obstacle coefficient to the discrete Laplacian eigenvalue
        # of the sampled sinusoidal equilibrium profile.  This tends to the
        # continuum pi^2/(2 eps^2) as dx/eps -> 0, while keeping an aligned
        # planar interface stationary at finite resolution.
        obstacle = 2.0 * np.sin(
            np.pi * cfg.grid_spacing / (2.0 * cfg.interface_width)
        ) ** 2 / cfg.grid_spacing**2
        # Sum_j Gamma [eta_j lap(eta_i) - eta_i lap(eta_j)
        #                    + pi^2/(2 eps^2) (eta_i - eta_j)].
        rate = cfg.intrinsic_mobility * cfg.gb_energy * (
            lap * phase_sum - eta_active * lap_sum
            + obstacle * (count * eta_active - phase_sum)
        ) * support
        if self.driving is not None:
            ext = np.asarray(self.driving(self.eta, self.time), dtype=float)
            if ext.shape != self.eta.shape:
                raise ValueError("driving callback returned the wrong shape")
            ext_active = ext[active_indices]
            ext_active -= (ext_active * support).sum(axis=0, keepdims=True) / count
            rate += cfg.intrinsic_mobility * ext_active * support
        rate *= self.mobility_scale[None, :, :]
        trial = eta_active + used_dt * rate
        # This is the bound/renormalization step used by the Qiu reference
        # implementation.  It preserves compact support unlike a smooth tail.
        np.clip(trial, 0.0, 1.0, out=trial)
        self.eta[active_indices] = trial / trial.sum(axis=0, keepdims=True)
        extinct = self.active_phases & (
            np.max(self.eta, axis=(1, 2)) < cfg.grain_extinction_threshold
        )
        if np.any(extinct) and np.count_nonzero(self.active_phases) > np.count_nonzero(extinct):
            self.active_phases[extinct] = False
            self.eta[extinct] = 0.0
            self.eta /= self.eta.sum(axis=0, keepdims=True)
        self.time += used_dt
        self.step_number += 1
        return StepDiagnostics(
            self.time, self.step_number, used_dt,
            free_energy(self.eta, cfg.gb_energy, cfg.interface_width, cfg.grid_spacing,
                        boundary=cfg.boundary_conditions),
            float(np.max(np.abs(self.eta.sum(axis=0) - 1.0))),
        )

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
                "active_phases": self.active_phases.copy()}

    def load_state_dict(self, state: dict[str, object]) -> None:
        self.eta = np.asarray(state["eta"], dtype=float).copy()
        self.time = float(state["time"])
        self.step_number = int(state["step_number"])
        self.mobility_scale = np.asarray(state.get("mobility_scale", np.ones(self.config.shape)), dtype=float).copy()
        self.active_phases = np.asarray(state.get("active_phases", np.max(self.eta, axis=(1, 2)) >= self.config.grain_extinction_threshold), dtype=bool).copy()
