from __future__ import annotations

from enum import Enum

import numpy as np

from grain_growth_pf.stochastic.hazard import CumulativeHazardClock


class ClimbStage(str, Enum):
    INACTIVE = "inactive"
    NUCLEATION = "nucleation"
    EXCHANGE = "exchange"
    TRANSPORT = "transport"
    COMPLETE = "quota_completion"


class SerialClimbCycle:
    """Explicit serial nucleation -> exchange -> transport state machine."""

    def __init__(self, rng: np.random.Generator, required_quota: float = 1.0):
        self.rng = rng
        self.required_quota = required_quota
        self.completed_quota = 0.0
        self.stage = ClimbStage.INACTIVE
        self.clock = CumulativeHazardClock(rng)
        self.history: list[tuple[float, ClimbStage]] = []
        self.last_completion_time: float | None = None
        self.last_transitions: list[tuple[float, ClimbStage, float]] = []

    def activate(self, time: float) -> None:
        if self.stage in {ClimbStage.INACTIVE, ClimbStage.COMPLETE}:
            self.stage = ClimbStage.NUCLEATION
            self.completed_quota = 0.0
            self.last_completion_time = None
            self.clock.reset()
            self.history.append((time, self.stage))

    def advance(self, dt: float, time: float, nucleation_rate: float,
                exchange_rate: float, transport_rate: float) -> bool:
        rates = {
            ClimbStage.NUCLEATION: nucleation_rate,
            ClimbStage.EXCHANGE: exchange_rate,
            ClimbStage.TRANSPORT: transport_rate,
        }
        self.last_transitions = []
        if dt < 0 or any(rate < 0 or not np.isfinite(rate) for rate in rates.values()):
            raise ValueError("serial rates and timestep must be finite and nonnegative")
        if self.stage not in rates:
            return self.stage == ClimbStage.COMPLETE
        remaining, current_time = float(dt), float(time)
        while remaining > 0 and self.stage in rates:
            rate = rates[self.stage]
            if rate == 0:
                self.clock.last_rate = 0.0
                break
            hazard_needed = max(float(self.clock.threshold) - self.clock.cumulative_hazard, 0.0)
            waiting_time = hazard_needed / rate
            if waiting_time > remaining:
                self.clock.cumulative_hazard += rate * remaining
                self.clock.last_rate = rate
                break
            self.clock.cumulative_hazard = float(self.clock.threshold)
            current_time += waiting_time
            remaining -= waiting_time
            if self.stage == ClimbStage.NUCLEATION:
                self.stage = ClimbStage.EXCHANGE
            elif self.stage == ClimbStage.EXCHANGE:
                self.stage = ClimbStage.TRANSPORT
            elif self.stage == ClimbStage.TRANSPORT:
                self.completed_quota = self.required_quota
                self.stage = ClimbStage.COMPLETE
                self.last_completion_time = current_time
            threshold = float(self.clock.threshold)
            self.history.append((current_time, self.stage))
            self.last_transitions.append((current_time, self.stage, threshold))
            self.clock.reset()
            if self.stage == ClimbStage.COMPLETE:
                return True
        return False

    @staticmethod
    def mean_completion_time(nucleation_rate: float, exchange_rate: float,
                             transport_rate: float) -> float:
        values = np.asarray((nucleation_rate, exchange_rate, transport_rate), float)
        return float("inf") if np.any(values <= 0) else float(np.sum(1.0 / values))
