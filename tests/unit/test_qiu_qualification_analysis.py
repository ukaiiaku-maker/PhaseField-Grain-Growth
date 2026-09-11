import json

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from analyze_qiu_full_field_qualification import (
    capture_assessment,
    compare_timestep_runs,
    load,
    pre_extinction_precursor,
    sha256,
    transition_timing,
)


def _trajectory(times: np.ndarray) -> pd.DataFrame:
    grain_size = 1.0 + times
    return pd.DataFrame({
        "time": times,
        "G_population": grain_size,
        "interfacial_energy": 10.0 + 2.0 * times,
        "elastic_energy": 0.5 + 0.1 * times,
        "compactness_mean": 1.1 + 0.01 * times,
        "compactness_p95": 1.2 + 0.01 * times,
        "compactness_max": 1.4 + 0.01 * times,
        "stress_p95": 0.02 + 0.002 * times,
        "grain_count": 100.0 - 5.0 * times,
        "largest_100_step_population_loss": np.zeros_like(times),
    })


def test_timestep_comparison_uses_time_progress_and_interpolated_histories():
    coarse = _trajectory(np.asarray([0.0, 1.0, 2.0]))
    refined = _trajectory(np.asarray([0.0, 0.5, 1.0, 1.5, 2.0]))

    result = compare_timestep_runs(coarse, refined)

    assert result["matched_time"]["time"] == 2.0
    assert result["matched_grain_size_progress"]["G_population"] == 3.0
    assert result["interpolated_energy_history"]["elastic_energy"]["normalized_linf"] < 1e-14
    assert result["all_gates_pass"]


def test_timestep_comparison_propagates_capture_classification():
    trajectory = _trajectory(np.asarray([0.0, 1.0, 2.0]))

    result = compare_timestep_runs(
        trajectory, trajectory.copy(), base_capture=True, fine_capture=False,
    )

    assert result["avalanche_indicators"]["corrected"]["detected"]
    assert not result["gates"]["same_avalanche_classification"]
    assert not result["all_gates_pass"]


def test_load_applies_external_source_attestation(tmp_path):
    diagnostic = tmp_path / "per_step_diagnostics.parquet"
    diagnostic.mkdir()
    pq.write_table(
        pa.Table.from_pydict({"step": [1], "time": [0.04]}),
        diagnostic / "part-000000.parquet",
    )
    (tmp_path / "manifest.json").write_text(json.dumps({"git_sha": "UNCOMMITTED"}))
    (tmp_path / "source_attestation.json").write_text(json.dumps({
        "verified_source_commit": "0123456789abcdef",
    }))

    _, manifest = load("corrected", tmp_path)

    assert manifest["git_sha"] == "UNCOMMITTED"
    assert manifest["verified_source_commit"] == "0123456789abcdef"


def test_capture_assessment_can_exclude_false_raw_guard_without_erasing_it(tmp_path):
    capture = tmp_path / "diagnostic_capture.json"
    capture.write_text(json.dumps({"step": 9001, "reason": "clipping_spike"}))
    (tmp_path / "diagnostic_capture_assessment.json").write_text(json.dumps({
        "raw_capture_sha256": sha256(capture),
        "accepted_as_transition": False,
        "scientific_transition_step": 9876,
        "reason": "Absolute clipping threshold fired on the ordinary obstacle halo.",
    }))

    result = capture_assessment(tmp_path)

    assert result["raw_capture_exists"]
    assert not result["accepted_as_transition"]
    assert result["scientific_transition_step"] == 9876
    assert result["raw_capture_sha256"] == sha256(capture)


def test_capture_assessment_requires_exact_raw_capture_hash(tmp_path):
    (tmp_path / "diagnostic_capture.json").write_text("{}")
    (tmp_path / "diagnostic_capture_assessment.json").write_text(json.dumps({
        "raw_capture_sha256": "0" * 64,
        "accepted_as_transition": False,
        "reason": "Known false marker.",
    }))

    with np.testing.assert_raises_regex(ValueError, "SHA-256 mismatch"):
        capture_assessment(tmp_path)


def test_transition_timing_separates_extinction_burst_and_morphology():
    step = np.arange(1, 8)
    frame = pd.DataFrame({
        "step": step,
        "grain_count": [100, 100, 99, 98, 88, 87, 87],
        "newly_extinct_phases": [0, 0, 1, 1, 10, 1, 0],
        "largest_100_step_population_loss": [0, 0, 1, 2, 12, 13, 13],
        "compactness_mean": [1.2, 1.2, 1.3, 1.4, 1.5, 2.6, 2.7],
        "compactness_max": [1.5, 1.5, 1.6, 1.7, 1.8, 2.0, 6.1],
        "disconnected_grain_count": [0, 0, 0, 1, 1, 2, 2],
        "total_energy": [10.0, 9.0, 8.0, 8.1, 7.0, 6.0, 5.0],
        "stress_linf": np.ones(7),
        "eigenstrain_linf": np.ones(7),
    })

    result = transition_timing(frame)

    assert result["first_population_decrease_step"] == 3
    assert result["first_extinction_step"] == 3
    assert result["first_disconnected_grain_step"] == 4
    assert result["first_100_step_population_loss_above_10pct_step"] == 5
    assert result["first_mean_compactness_above_2p5_step"] == 6
    assert result["first_max_compactness_above_6_step"] == 7
    assert result["first_complete_energy_increase_step"] == 4
    assert result["first_objective_avalanche_indicator_step"] == 5


def test_pre_extinction_precursor_uses_only_rows_before_first_extinction():
    step = np.arange(1, 7)
    frame = pd.DataFrame({
        "step": step,
        "grain_count": [10, 10, 10, 10, 9, 8],
        "newly_extinct_phases": [0, 0, 0, 0, 1, 1],
        "disconnected_grain_count": np.zeros(6),
        "compactness_max": np.ones(6),
        "source_increment_l2": [1, 2, 3, 4, 100, 100],
        "stress_linf": [1, 2, 3, 4, 100, 100],
        "eigenstrain_linf": [1, 2, 3, 4, 100, 100],
        "interfacial_energy": [1, 2, 3, 4, 100, 100],
        "elastic_energy": [1, 2, 3, 4, 100, 100],
        "total_energy": [1, 2, 3, 4, 100, 100],
    })

    result = pre_extinction_precursor(frame)

    assert result["first_extinction_step"] == 5
    assert result["step_last"] == 4
    assert result["source_increment_l2_median_change"] == 2.0
