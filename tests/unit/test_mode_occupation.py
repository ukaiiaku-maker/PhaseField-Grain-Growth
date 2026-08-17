import json

import pandas as pd
import pytest

from grain_growth_pf.analysis.mode_occupation import analyze_mode_occupation


def _write_run(root, regime, temperature, seed, rows):
    run = root / f"{regime}-T{temperature:g}-s{seed}"
    run.mkdir()
    (run / "manifest.json").write_text(json.dumps({
        "status": "completed",
        "git_sha": "simulation-sha",
        "config": {
            "regime": regime,
            "seed": seed,
            "pf": {"temperature": temperature},
        },
    }))
    pd.DataFrame(rows).to_csv(run / "events.csv", index=False)
    return run


def test_mode_occupation_uses_completed_modes_and_bootstraps_runs(tmp_path):
    rows = {
        "event_type": [
            "activation_hit", "disconnection_mode", "disconnection_mode",
            "compatibility_release",
        ],
        "barrier_type": ["easy", "easy", "intermediate", "high"],
        "activation_volume": [None, "0.25;0.25", "-0.25;0.5", "0.25;1.0"],
        "shear_strain_increment": [0.0, 0.1, 0.2, 0.4],
        "volumetric_strain_increment": [0.0, 0.01, 0.02, 0.04],
    }
    runs = [
        _write_run(tmp_path, "SC3", temperature, seed, rows)
        for temperature in (800.0, 900.0)
        for seed in (1, 2)
    ]
    (tmp_path / "campaign_manifest.json").write_text(json.dumps({
        "status": "completed", "runs": [str(run) for run in runs]
    }))

    result = analyze_mode_occupation(tmp_path, "SC3", bootstrap_samples=20)

    assert result["simulation_git_shas"] == ["simulation-sha"]
    assert len(result["temperature_summaries"]) == 2
    summary = result["temperature_summaries"][0]
    assert summary["mode_events"] == 6
    assert summary["mean_signed_beta"] == 1.0
    assert summary["mean_absolute_beta"] == 7 / 3
    assert summary["absolute_beta_fractions"] == {
        "1": 1 / 3, "2": 1 / 3, "4": 1 / 3,
    }
    assert summary["family_fractions"] == {
        "easy": 1 / 3, "high": 1 / 3, "intermediate": 1 / 3,
    }
    assert summary["signed_shear_strain"] == pytest.approx(1.4)
    assert summary["mean_absolute_beta_95pct"] == [7 / 3, 7 / 3]
