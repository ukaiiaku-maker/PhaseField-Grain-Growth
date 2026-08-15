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
        if self.stage not in rates:
            return self.stage == ClimbStage.COMPLETE
        events = self.clock.advance(rates[self.stage], dt, time)
        for event in events:
            if self.stage == ClimbStage.NUCLEATION:
                self.stage = ClimbStage.EXCHANGE
            elif self.stage == ClimbStage.EXCHANGE:
                self.stage = ClimbStage.TRANSPORT
            elif self.stage == ClimbStage.TRANSPORT:
                self.completed_quota = self.required_quota
                self.stage = ClimbStage.COMPLETE
                self.last_completion_time = event.event_time
            self.history.append((event.event_time, self.stage))
            self.clock.reset()
            if self.stage == ClimbStage.COMPLETE:
                return True
        return False

    @staticmethod
    def mean_completion_time(nucleation_rate: float, exchange_rate: float,
                             transport_rate: float) -> float:
        values = np.asarray((nucleation_rate, exchange_rate, transport_rate), float)
        return float("inf") if np.any(values <= 0) else float(np.sum(1.0 / values))
