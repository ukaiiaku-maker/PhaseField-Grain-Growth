from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from numpy.typing import NDArray

from grain_growth_pf.config import PFConfig
from .free_energy import free_energy
from .kernels import pairwise_obstacle_step, pairwise_obstacle_trial_statistics

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
    requested_dt: float = float("nan")
    max_abs_order_parameter_increment: float = float("nan")
    max_raw_order_parameter_increment: float = float("nan")
    clipped_low_count: int = 0
    clipped_high_count: int = 0
    trial_value_count: int = 0
    clipped_fraction: float = 0.0
    total_renormalization_correction: float = 0.0
    newly_extinct_phases: int = 0
    minimum_pre_extinction_maximum: float = float("nan")
    max_external_pair_driving_difference: float = 0.0


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
        sums = eta.sum(axis=0)
        already_on_simplex = bool(
            np.all(np.isfinite(eta))
            and np.min(eta) >= 0.0
            and np.max(np.abs(sums - 1.0)) <= 1e-12
        )
        self.eta = eta.copy() if already_on_simplex else project_simplex(eta)
        self.active_phases = np.max(self.eta, axis=(1, 2)) >= config.grain_extinction_threshold
        if not np.any(self.active_phases):
            raise ValueError("initial condition contains no active grain")
        self.config = config
        self.driving = driving
        self.mobility_scale = np.ones(config.shape, dtype=float)
        self.time = 0.0
        self.step_number = 0
        self.capture_step_diagnostics = False

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
        before_eta = self.eta
        before_active = self.active_phases.copy() if self.capture_step_diagnostics else self.active_phases
        self.eta = pairwise_obstacle_step(
            before_eta,
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
        raw = None
        if self.capture_step_diagnostics:
            raw = pairwise_obstacle_trial_statistics(
                before_eta, before_active, self.mobility_scale, external,
                use_external, used_dt, cfg.intrinsic_mobility, cfg.gb_energy,
                cfg.interface_width, cfg.grid_spacing,
                cfg.boundary_conditions == "periodic",
            )
        extinct = self.active_phases & (
            np.max(self.eta, axis=(1, 2)) < cfg.grain_extinction_threshold
        )
        pre_extinction_maximum = (
            float(np.min(np.max(before_eta[extinct], axis=(1, 2))))
            if self.capture_step_diagnostics and np.any(extinct) else float("nan")
        )
        if np.any(extinct) and np.count_nonzero(self.active_phases) > np.count_nonzero(extinct):
            self.active_phases[extinct] = False
            self.eta[extinct] = 0.0
            self.eta /= self.eta.sum(axis=0, keepdims=True)
        self.time += used_dt
        self.step_number += 1
        return StepDiagnostics(
            self.time, self.step_number, used_dt,
            (
                free_energy(
                    self.eta, cfg.gb_energy, cfg.interface_width, cfg.grid_spacing,
                    boundary=cfg.boundary_conditions,
                )
                if compute_energy else float("nan")
            ),
            float(np.max(np.abs(self.eta.sum(axis=0) - 1.0))),
            requested_dt=float(requested),
            max_abs_order_parameter_increment=(
                float(np.max(np.abs(self.eta - before_eta)))
                if self.capture_step_diagnostics else float("nan")
            ),
            max_raw_order_parameter_increment=(
                float(raw["max_raw_increment"]) if raw is not None else float("nan")
            ),
            clipped_low_count=int(raw["clipped_low"]) if raw is not None else 0,
            clipped_high_count=int(raw["clipped_high"]) if raw is not None else 0,
            trial_value_count=int(raw["trial_count"]) if raw is not None else 0,
            clipped_fraction=(
                float((raw["clipped_low"] + raw["clipped_high"]) / raw["trial_count"])
                if raw is not None and raw["trial_count"] else 0.0
            ),
            total_renormalization_correction=(
                float(raw["renormalization_l1"]) if raw is not None else 0.0
            ),
            newly_extinct_phases=int(np.count_nonzero(extinct)),
            minimum_pre_extinction_maximum=pre_extinction_maximum,
            max_external_pair_driving_difference=(
                float(raw["max_external_pair_difference"]) if raw is not None else 0.0
            ),
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
