from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


_SPEC = importlib.util.spec_from_file_location(
    "analyze_shear_screen",
    Path(__file__).parents[2] / "scripts" / "analyze_shear_screen.py",
)
assert _SPEC is not None and _SPEC.loader is not None
_ANALYSIS = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_ANALYSIS)


def test_topology_window_metrics_use_matched_population_not_equal_time():
    growth = pd.DataFrame({
        "step": [0, 10, 20, 30, 40],
        "time": [0.0, 1.0, 2.0, 3.0, 4.0],
        "mean_radius": [1.0, 1.1, 1.2, 1.3, 1.4],
        "grain_count": [200, 190, 180, 170, 160],
        "radius_rate": [0.1] * 5,
    })
    events = pd.DataFrame({
        "step": [12, 22, 32],
        "event_type": [
            "compatibility_release", "tj_compatibility_release", "gb_sink_completion",
        ],
    })
    first = _ANALYSIS._window_metrics("case", growth, events)[0]
    assert first["available"]
    assert np.isclose(first["Rdot"], 0.1)
    assert np.isclose(first["f_G"], 1 / 3)
    assert np.isclose(first["f_T"], 1 / 3)
    assert np.isclose(first["f_C"], 1 / 3)


def test_solver_step_response_uses_exact_entity_steps_and_all_windows():
    steps = np.arange(65)
    growth = pd.DataFrame({
        "step": steps,
        "time": steps * 0.04,
        "mean_radius": 1.0 + steps * 0.001,
        "grain_count": np.linspace(190, 160, len(steps)),
        "radius_rate": np.full(len(steps), 0.025),
    })
    trace = pd.DataFrame({
        "step": steps,
        "entity_id": "gb:1-2:0",
        "local_normal_velocity": steps.astype(float),
        "shear_state_s": 64.0 - steps,
        "tau_int": 2.0 * (64.0 - steps),
    })
    index = pd.DataFrame({
        "event_id": ["event-1"],
        "event_type": ["compatibility_release"],
        "sink_path": [""],
        "event_step": [32],
        "entity_id": ["gb:1-2:0"],
    })
    detail, nulls = _ANALYSIS._event_responses("case", growth, trace, index)
    assert {row["window_steps"] for row in detail} == {1, 2, 4, 8, 16, 32}
    assert {row["window_steps"] for row in nulls} == {1, 2, 4, 8, 16, 32}
    one_step = next(row for row in detail if row["window_steps"] == 1)
    assert one_step["signed_velocity_response"] == 2.0
    assert one_step["shear_state_relaxation"] == -2.0
