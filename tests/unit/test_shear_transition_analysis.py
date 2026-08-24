from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np


_SPEC = importlib.util.spec_from_file_location(
    "analyze_shear_transition",
    Path(__file__).parents[2] / "scripts" / "analyze_shear_transition.py",
)
assert _SPEC is not None and _SPEC.loader is not None
_ANALYSIS = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_ANALYSIS)


def test_weighted_quantile_respects_length_weights():
    value = _ANALYSIS.weighted_quantile(
        np.asarray([0.0, 1.0, 2.0]), np.asarray([1.0, 8.0, 1.0]), 0.5
    )
    assert np.isclose(value, 1.0)


def test_arrest_episodes_break_on_trace_gaps_and_report_curvature():
    rows = _ANALYSIS.arrest_episodes(
        np.asarray([1, 2, 3, 7, 8, 9]),
        np.asarray([False, True, True, True, True, False]),
        np.asarray([0.1, 0.2, 0.3, 0.7, 0.8, 0.9]),
        0.04,
    )
    assert [row["duration_steps"] for row in rows] == [2, 2]
    assert np.isclose(rows[0]["entry_curvature"], 0.2)
    assert np.isclose(rows[1]["exit_curvature"], 0.8)


def test_segmented_change_point_localizes_known_hinge():
    x = np.asarray([0.10, 0.20, 0.25, 0.30, 0.35, 0.40])
    y = 1.0 - x - 8.0 * np.maximum(0.0, x - 0.30)
    fit = _ANALYSIS.segmented_change_point(x, y)
    assert np.isclose(fit["critical_Ks"], 0.30)
    assert fit["sse"] < 1e-20
