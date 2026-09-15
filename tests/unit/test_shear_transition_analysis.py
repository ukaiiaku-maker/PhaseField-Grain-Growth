from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


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


def test_frame_windows_separates_exact_domain_and_length_weighted_force_metrics(tmp_path):
    run = tmp_path / "case"
    frames = run / "frames"
    frames.mkdir(parents=True)
    np.savez_compressed(
        frames / "frame-0000000.npz",
        grain_count=np.asarray(180),
        boundary_mask=np.asarray([[1, 1], [1, 0]], dtype=np.uint8),
        shear=np.asarray([[1.0, 1.0], [1.0, 0.0]]),
        shear_stress=np.asarray([[-0.2, -0.3], [-0.4, 0.0]]),
        p_shear=np.asarray([[-0.1, -0.2], [-0.3, 0.0]]),
        chi_s=np.asarray([[0.5, 0.9], [1.1, np.nan]]),
        p_net=np.asarray([[1.0, -0.2], [-0.4, np.nan]]),
        active_shear_length=np.asarray(3.0),
        total_gb_length=np.asarray(3.0),
        stored_shear_energy=np.asarray(1.5),
        stored_shear_energy_per_active_gb_length=np.asarray(0.5),
        mean_abs_tau_int_active_domain=np.asarray(7.0),
        fraction_active_gb_length_chi_s_near_one=np.asarray(0.625),
        fraction_active_gb_length_chi_s_gt_one=np.asarray(0.25),
    )
    rows = _ANALYSIS._frame_windows("GTSC_GB_Ks025", 5101, 0.25, run, 1.0)
    row = next(item for item in rows if item["topology_window"] == "N190_to_160")
    assert row["available"]
    assert np.isclose(row["mean_abs_tau_int_active_domain"], 7.0)
    assert np.isclose(row["mean_abs_tau_int_active_length"], 0.3)
    assert np.isclose(row["chi_s_p50_active_length"], 0.9)
    assert np.isclose(row["fraction_active_gb_length_chi_s_near_one"], 0.625)
    assert np.isclose(row["fraction_active_gb_length_chi_s_gt_one"], 0.25)
    assert np.isclose(row["fraction_active_gb_length_p_net_negative"], 2.0 / 3.0)


def test_experimental_observables_reports_bursts_waits_and_grain_velocity(tmp_path):
    growth = pd.DataFrame({
        "step": [0, 1, 2, 3], "time": [0.0, 1.0, 2.0, 3.0],
        "grain_count": [190, 180, 170, 160], "mean_radius": [1.0, 1.0, 2.0, 2.0],
    })
    run = tmp_path / "case"
    run.mkdir()
    pd.DataFrame({
        "time": [0.0, 1.0, 2.0, 3.0], "step": [0, 1, 2, 3],
        "grain_id": [7, 7, 7, 7], "radius": [1.0, 1.5, 2.0, 2.5],
    }).to_csv(run / "grain_tracks.csv", index=False)
    rows = _ANALYSIS._experimental_observables(
        "GTSC_GB_Ks025", 5101, 0.25, growth, run,
    )
    row = next(item for item in rows if item["topology_window"] == "N190_to_160")
    assert row["available"]
    assert np.isclose(row["stationary_fraction"], 2.0 / 3.0)
    assert np.isclose(row["fraction_positive_growth_in_largest_10pct_bursts"], 1.0)
    assert row["waiting_episodes"] == 2
    assert np.isclose(row["waiting_time_median"], 1.0)
    assert row["grain_velocity_samples"] == 3
    assert np.isclose(row["grain_abs_radial_velocity_p95"], 0.5)
