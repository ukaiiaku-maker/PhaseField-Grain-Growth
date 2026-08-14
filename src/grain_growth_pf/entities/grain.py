from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Grain:
    grain_id: int
    orientation: float
    area: float = 0.0
    equivalent_radius: float = 0.0
    centroid: tuple[float, float] = (0.0, 0.0)
    neighbors: set[int] = field(default_factory=set)
    perimeter: float = 0.0
    event_history: list[str] = field(default_factory=list)

