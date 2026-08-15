from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import gammainc

from .hazard import CumulativeHazardClock, HazardEvent


def poisson_completion_probability(required_hits: int, window_hazard: float) -> float:
    if required_hits < 1 or window_hazard < 0:
        raise ValueError("K >= 1 and Lambda >= 0 are required")
    # P[N >= K] = gammainc(K, Lambda); stable at both tails.
    return float(gammainc(required_hits, window_hazard))


@dataclass
class CompletionEvent:
    time: float
    hits: int


class MultiHitProcess:
    def __init__(self, required_hits: int, rng: np.random.Generator,
                 interpretation: str = "persistent_hits"):
        if required_hits < 1:
            raise ValueError("required_hits must be at least one")
        if interpretation not in {"persistent_hits", "packet_reset"}:
            raise ValueError("unknown multihit interpretation")
        self.required_hits = required_hits
        self.interpretation = interpretation
        self.hit_count = 0
        self.clock = CumulativeHazardClock(rng)
        self.last_hit_events: list[HazardEvent] = []
        self.last_hit_counts: list[int] = []
        self.last_hit_completions: list[bool] = []

    def begin_window(self) -> None:
        """Start an encounter packet while applying its declared memory rule."""
        if self.interpretation == "packet_reset":
            self.hit_count = 0
        self.clock.reset()

    def advance(self, rate: float, dt: float, time: float) -> list[CompletionEvent]:
        completions: list[CompletionEvent] = []
        self.last_hit_events = self.clock.advance(rate, dt, time)
        self.last_hit_counts = []
        self.last_hit_completions = []
        for event in self.last_hit_events:
            self.hit_count += 1
            completed = self.hit_count >= self.required_hits
            self.last_hit_counts.append(self.hit_count)
            self.last_hit_completions.append(completed)
            if completed:
                completions.append(CompletionEvent(event.event_time, self.hit_count))
                self.hit_count = 0
        return completions

    def close_window(self) -> bool:
        complete = self.hit_count >= self.required_hits
        if self.interpretation == "packet_reset":
            self.hit_count = 0
            self.clock.reset()
        return complete
