from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from numpy.typing import NDArray

from grain_growth_pf.config import PFConfig
from .free_energy import chemical_potential, free_energy

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
    """Constrained Allen-Cahn multiphase-field reference solver.

    The Lagrange multiplier is the local mean chemical potential, ensuring
    sum(d eta_i/dt)=0 before the bound-preserving simplex projection. External
    pair physics enters as a zero-sum phase driving field.
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
        mu = chemical_potential(
            self.eta, cfg.gb_energy, cfg.interface_width,
            cfg.grid_spacing, cfg.boundary_conditions,
        )
        # A globally active phase may advance by one stencil cell beyond its
        # current diffuse support. Restricting the rate to eta>0 exactly pins
        # neighbor switches because a phase can never enter a zero-valued
        # adjacent pixel. The one-cell dilation permits local topology changes
        # without allowing remote grain nucleation.
        present = self.eta > 1e-12
        active = present.copy()
        for axis in (-2, -1):
            forward = np.roll(present, 1, axis=axis)
            backward = np.roll(present, -1, axis=axis)
            if cfg.boundary_conditions != "periodic":
                first = [slice(None)] * present.ndim
                last = [slice(None)] * present.ndim
                first[axis] = 0
                last[axis] = -1
                forward[tuple(first)] = False
                backward[tuple(last)] = False
            active |= forward | backward
        active &= self.active_phases[:, None, None]
        count = np.maximum(active.sum(axis=0, keepdims=True), 1)
        lagrange = (mu * active).sum(axis=0, keepdims=True) / count
        # For the chosen equilibrium profile integral(|grad eta|^2) = 1/(3w).
        # L=M_sharp/(3w) therefore gives v_n=M_sharp*gamma*kappa.
        kinetic = cfg.intrinsic_mobility / (3.0 * cfg.interface_width)
        rate = -kinetic * (mu - lagrange) * active
        if self.driving is not None:
            ext = np.asarray(self.driving(self.eta, self.time), dtype=float)
            if ext.shape != self.eta.shape:
                raise ValueError("driving callback returned the wrong shape")
            ext -= ext.mean(axis=0, keepdims=True)
            rate += kinetic * ext
        rate *= self.mobility_scale[None, :, :]
        trial = self.eta + used_dt * rate
        if float(trial.min()) >= 0.0 and float(trial.max()) <= 1.0:
            # The constrained rate has zero local sum. Normalization removes
            # only roundoff and avoids an O(N log N) simplex sort at every pixel.
            self.eta = trial / trial.sum(axis=0, keepdims=True)
        else:
            self.eta = project_simplex(trial)
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
            free_energy(self.eta, cfg.gb_energy, cfg.interface_width, cfg.grid_spacing),
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
