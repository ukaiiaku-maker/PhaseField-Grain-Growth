from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

K_B_EV = 8.617333262145e-5


@dataclass(frozen=True)
class ModeDriving:
    normal_pressure: float = 0.0
    resolved_shear: float = 0.0
    vacancy_chemical_potential: float = 0.0


@dataclass(frozen=True)
class DisconnectionMode:
    mode_id: str
    burgers: tuple[float, float]
    step_height: float
    point_defect_quota: float
    barrier_ev: float
    attempt_frequency: float
    site_multiplicity: float
    activation_volume_normal: float = 0.0
    activation_volume_shear: float = 0.0
    activation_vacancies: float = 0.0
    delta_s: float = 0.0
    delta_q: float = 0.0
    family: str = "easy"

    def __post_init__(self) -> None:
        if self.attempt_frequency < 0 or self.site_multiplicity < 0 or self.barrier_ev < 0:
            raise ValueError("frequency, multiplicity, and zero-driving barrier must be nonnegative")
        if self.step_height == 0:
            raise ValueError("zero step height is not a migrating disconnection mode")

    @property
    def beta(self) -> float:
        tangent_b = float(np.linalg.norm(self.burgers))
        return tangent_b / self.step_height

    def resolved_shear(self, stress: NDArray[np.float64], normal: NDArray[np.float64]) -> float:
        n = np.asarray(normal, dtype=float)
        n /= max(np.linalg.norm(n), np.finfo(float).tiny)
        b = np.asarray(self.burgers, dtype=float)
        if np.linalg.norm(b) == 0:
            return 0.0
        t = b / np.linalg.norm(b)
        return float(t @ np.asarray(stress, dtype=float) @ n)

    def activation_work_ev(self, driving: ModeDriving) -> float:
        return (
            driving.normal_pressure * self.activation_volume_normal
            + driving.resolved_shear * self.activation_volume_shear
            + driving.vacancy_chemical_potential * self.activation_vacancies
        )

    def effective_barrier_ev(self, driving: ModeDriving) -> float:
        # Defined high-driving limit: once work reaches G0, the transition
        # state is barrierless and the process becomes attempt limited.
        return max(0.0, self.barrier_ev - self.activation_work_ev(driving))

    def log_rate(self, temperature: float, driving: ModeDriving = ModeDriving()) -> float:
        if temperature <= 0:
            raise ValueError("temperature must be positive")
        prefactor = self.site_multiplicity * self.attempt_frequency
        if prefactor == 0:
            return -np.inf
        return float(np.log(prefactor) - self.effective_barrier_ev(driving) / (K_B_EV * temperature))

    def rate(self, temperature: float, driving: ModeDriving = ModeDriving()) -> float:
        log_r = self.log_rate(temperature, driving)
        max_log = np.log(self.site_multiplicity * self.attempt_frequency) if self.site_multiplicity * self.attempt_frequency else -np.inf
        return float(np.exp(min(log_r, max_log)))

