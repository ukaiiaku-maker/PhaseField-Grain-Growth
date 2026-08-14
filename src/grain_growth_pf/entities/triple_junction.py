from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class TripleJunction:
    grain_ids: tuple[int, int, int]
    position: tuple[float, float]
    travel_distance: float = 0.0
    adjoining_boundaries: set[str] = field(default_factory=set)
    compatible: bool = True
    residual_burgers: np.ndarray = field(default_factory=lambda: np.zeros(2))
    event_history: list[str] = field(default_factory=list)
    age: int = 0

    @property
    def entity_id(self) -> str:
        return "tj:" + "-".join(map(str, sorted(self.grain_ids)))

    def add_burgers(self, increment: np.ndarray) -> None:
        value = np.asarray(increment, dtype=float)
        if value.shape != (2,):
            raise ValueError("TJ Burgers increment must be a 2-vector")
        self.residual_burgers += value

