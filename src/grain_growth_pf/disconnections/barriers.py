from __future__ import annotations

from dataclasses import replace
from typing import Iterable

import numpy as np

from .mode import DisconnectionMode


def assign_barriers(modes: Iterable[DisconnectionMode], kind: str, seed: int,
                    mean_ev: float, std_ev: float = 0.0,
                    bounds_ev: tuple[float, float] | None = None,
                    misorientation: float | None = None,
                    character_coefficient_ev: float = 0.0) -> list[DisconnectionMode]:
    """Return a quenched barrier realization attached to copied mode objects."""
    rng = np.random.default_rng(seed)
    source = list(modes)
    if kind == "fixed":
        values = np.full(len(source), mean_ev)
    elif kind == "truncated_gaussian":
        if bounds_ev is None:
            raise ValueError("truncated_gaussian requires bounds_ev")
        values = rng.normal(mean_ev, std_ev, len(source))
        values = np.clip(values, *bounds_ev)
    elif kind == "bounded_two_level":
        if bounds_ev is None:
            raise ValueError("bounded_two_level requires bounds_ev")
        values = rng.choice(np.asarray(bounds_ev), len(source))
    elif kind == "gb_character":
        if misorientation is None:
            raise ValueError("gb_character requires misorientation")
        values = np.full(len(source), mean_ev + character_coefficient_ev * abs(np.sin(misorientation)))
    else:
        raise ValueError(f"unknown barrier distribution {kind}")
    return [replace(mode, barrier_ev=float(max(value, 0.0))) for mode, value in zip(source, values)]


def renew_barrier(mode: DisconnectionMode, rng: np.random.Generator, std_ev: float,
                  bounds_ev: tuple[float, float]) -> DisconnectionMode:
    """Explicit event-to-event (annealed) barrier renewal."""
    value = float(np.clip(rng.normal(mode.barrier_ev, std_ev), *bounds_ev))
    return replace(mode, barrier_ev=value)

