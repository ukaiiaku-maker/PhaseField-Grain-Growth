from __future__ import annotations

import numpy as np
import pandas as pd

from grain_growth_pf.analysis.mechanism_diagnostics import (
    event_burst_coupling,
    factorial_effects,
    fit_pin_lifetime_models,
    pin_episodes,
    reconstruct_gb_occupancy,
)


def test_event_burst_coupling_exceeds_grain_stratified_shuffle() -> None:
    rows = []
    events = []
    for grain_id in range(12):
        radius = 1.0
        for time in range(21):
            if time and time in (5, 10, 15, 20):
                radius += 1.0
            rows.append({"grain_id": grain_id, "time": float(time), "radius": radius})
        for time in (5, 10, 15, 20):
            events.append({
                "time": float(time), "event_type": "compatibility_release",
                "grain_ids": str(grain_id), "entity_id": f"gb:{grain_id}-99:0",
            })
    result = event_burst_coupling(
        pd.DataFrame(rows), pd.DataFrame(events), shuffle_samples=64, seed=8
    )
    assert result["large_burst_risk_ratio_excess"] > 1.5
    assert result["top_10pct_release_excess"] > 0.5


def test_event_free_control_has_zero_excess() -> None:
    tracks = pd.DataFrame({
        "grain_id": [1, 1, 1], "time": [0.0, 1.0, 2.0],
        "radius": [1.0, 1.0, 2.0],
    })
    result = event_burst_coupling(tracks, pd.DataFrame(), shuffle_samples=8)
    assert result["release_count"] == 0
    assert result["top_5pct_release_excess"] == 0.0


def test_factorial_effects_recovers_main_and_interaction() -> None:
    rows = []
    factors = ("G", "T", "S", "C")
    for mask in range(16):
        active = {factor for index, factor in enumerate(factors) if mask & (1 << index)}
        regime = "B0" if not active else "".join(factor for factor in factors if factor in active)
        x = {factor: (1.0 if factor in active else -1.0) for factor in factors}
        response = 7.0 + 2.5 * x["C"] - 1.25 * x["S"] * x["C"]
        rows.append({"regime": regime, "response": response})
    effects = factorial_effects(pd.DataFrame(rows), ["response"]).set_index("effect")
    assert np.isclose(effects.loc["C", "value"], 5.0)
    assert np.isclose(effects.loc["SxC", "value"], -2.5)
    assert np.isclose(effects.loc["G", "value"], 0.0)


def test_pin_survival_and_gc_occupancy_reconstruction() -> None:
    boundaries = pd.DataFrame({
        "entity_id": ["gb:1-2:0"] * 6,
        "time": np.arange(6.0),
        "blocked": [False, True, True, True, False, False],
    })
    episodes = pin_episodes(boundaries)
    assert len(episodes) == 1
    assert episodes.iloc[0]["duration"] == 3.0
    events = pd.DataFrame({
        "entity_id": ["gb:1-2:0", "gb:1-2:0"],
        "time": [2.0, 4.0],
        "event_type": ["compatibility_release", "climb_quota_completion"],
    })
    occupancy = reconstruct_gb_occupancy(episodes, events, has_g=True, has_c=True)
    assert np.isclose(occupancy["multiple"], 1.0)
    assert np.isclose(occupancy["C_limited"], 2.0)

    c_only_events = events.iloc[[1]]
    c_only = reconstruct_gb_occupancy(
        episodes, c_only_events, has_g=True, has_c=True
    )
    assert np.isclose(c_only["C_limited"], 3.0)
    assert c_only["multiple"] == 0.0


def test_censored_lifetime_fits_are_finite() -> None:
    episodes = pd.DataFrame({
        "duration": [1.0, 2.0, 3.0, 4.0, 6.0],
        "censored": [False, False, False, True, False],
    })
    fitted = fit_pin_lifetime_models(episodes)
    assert set(fitted["model"]) == {"exponential", "gamma", "weibull"}
    assert np.all(np.isfinite(fitted["AIC"]))
