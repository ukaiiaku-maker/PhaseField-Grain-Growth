import json

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from analyze_qiu_full_field_qualification import compare_timestep_runs, load


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
