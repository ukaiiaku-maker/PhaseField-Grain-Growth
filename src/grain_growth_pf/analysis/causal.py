from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np


def _window_delta(profile: Mapping[int, float], center: int) -> float | None:
    before = [profile.get(center + offset) for offset in range(-5, 0)]
    after = [profile.get(center + offset) for offset in range(1, 6)]
    if any(value is None or not np.isfinite(value) for value in (*before, *after)):
        return None
    return float(np.mean(after) - np.mean(before))


def trace_center_causal_null(
    profiles: Sequence[Mapping[int, float]],
    *,
    shuffles: int = 200,
    seed: int = 5101,
) -> dict[str, float | int]:
    """Compare release-centered motion with within-trace pseudo-event centers.

    Each profile is already conditioned on a deterministically sampled event.
    The null retains that event's local entities and velocity history, but moves
    the alignment center away from the true event.  It therefore tests temporal
    alignment, not whether event-bearing neighborhoods differ from the bulk.
    """
    usable: list[tuple[float, list[float]]] = []
    for profile in profiles:
        actual = _window_delta(profile, 0)
        if actual is None:
            continue
        centers = [
            center
            for center in range(min(profile) + 5, max(profile) - 4)
            if abs(center) > 5 and _window_delta(profile, center) is not None
        ]
        if not centers:
            continue
        usable.append((actual, [_window_delta(profile, center) for center in centers]))
    if not usable:
        return {
            "events": 0,
            "actual_mean_delta": float("nan"),
            "null_mean_delta": float("nan"),
            "causal_excess": float("nan"),
            "null_sd": float("nan"),
            "p_one_sided": float("nan"),
            "shuffles": shuffles,
        }
    rng = np.random.default_rng(seed)
    actual_mean = float(np.mean([item[0] for item in usable]))
    null = np.empty(shuffles)
    for index in range(shuffles):
        null[index] = np.mean(
            [choices[int(rng.integers(0, len(choices)))] for _, choices in usable]
        )
    return {
        "events": len(usable),
        "actual_mean_delta": actual_mean,
        "null_mean_delta": float(null.mean()),
        "causal_excess": float(actual_mean - null.mean()),
        "null_sd": float(null.std(ddof=1)),
        "p_one_sided": float((1 + np.count_nonzero(null >= actual_mean)) / (shuffles + 1)),
        "shuffles": shuffles,
    }
