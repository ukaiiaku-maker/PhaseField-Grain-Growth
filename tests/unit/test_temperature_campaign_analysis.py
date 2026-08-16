import json

import numpy as np
import pandas as pd
import pytest

from grain_growth_pf.analysis.campaign import analyze_campaign
from grain_growth_pf.config import ModelConfig, PFConfig
from grain_growth_pf.disconnections.mode import K_B_EV


def test_temperature_campaign_uses_one_common_exponent(tmp_path):
    temperatures = np.array([750.0, 850.0, 950.0, 1100.0])
    activation_energy = 0.42
    coefficients = 100.0 * np.exp(-activation_energy / (K_B_EV * temperatures))
    run_paths = []
    for temperature, coefficient in zip(temperatures, coefficients):
        for seed in (1, 2):
            run_dir = tmp_path / f"T{temperature:g}-s{seed}"
            run_dir.mkdir()
            config = ModelConfig(
                regime="synthetic", seed=seed,
                pf=PFConfig(shape=(8, 8), temperature=temperature),
            )
            (run_dir / "manifest.json").write_text(json.dumps({
                "config": config.to_dict(), "git_sha": "synthetic-sha", "status": "completed"
            }))
            time = np.linspace(0.0, 20.0, 31)
            radius = np.sqrt(16.0 + coefficient * time)
            pd.DataFrame({
                "run_id": f"s{seed}", "time": time, "step": np.arange(len(time)),
                "grain_id": np.ones(len(time), dtype=int), "area": np.pi * radius**2,
                "radius": radius, "perimeter": 2.0 * np.pi * radius,
                "neighbors": np.full(len(time), 6),
            }).to_csv(run_dir / "grain_tracks.csv", index=False)
            run_paths.append(str(run_dir))
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    (campaign / "campaign_manifest.json").write_text(json.dumps({
        "runs": run_paths, "status": "completed",
    }))

    summary = analyze_campaign(campaign, bootstrap_samples=5)
    assert np.allclose(summary["n"], 2.0, atol=1e-6)
    assert np.allclose(summary.sort_values("temperature")["K"], coefficients, rtol=1e-6)
    assert np.allclose(summary["Q_app"], activation_energy, atol=1e-6)
    diagnostics = json.loads((campaign / "mechanism_summary_diagnostics.json").read_text())
    assert all("temperature_series_fit" in detail for detail in diagnostics)
    assert all(detail["provenance"]["campaign"] == str(campaign) for detail in diagnostics)
    assert all(
        detail["provenance"]["simulation_git_shas"] == ["synthetic-sha"]
        for detail in diagnostics
    )
    (campaign / "campaign_manifest.json").write_text(json.dumps({
        "runs": run_paths, "status": "running",
    }))
    with pytest.raises(ValueError, match="campaign is not complete"):
        analyze_campaign(campaign, bootstrap_samples=1)
    diagnostic_runs = run_paths + [str(campaign / "not-started")]
    (campaign / "campaign_manifest.json").write_text(json.dumps({
        "runs": diagnostic_runs, "status": "running",
    }))
    diagnostic = analyze_campaign(
        campaign, bootstrap_samples=1, require_completed=False
    )
    assert len(diagnostic) == len(temperatures)


def test_stagnant_temperature_series_suppresses_common_exponent_and_activation(tmp_path):
    run_paths = []
    for temperature_index, temperature in enumerate((800.0, 900.0, 1000.0, 1100.0)):
        for seed in (1, 2):
            run_dir = tmp_path / f"flat-T{temperature:g}-s{seed}"
            run_dir.mkdir()
            config = ModelConfig(
                regime="flat", seed=seed,
                pf=PFConfig(shape=(8, 8), temperature=temperature),
            )
            (run_dir / "manifest.json").write_text(json.dumps({
                "config": config.to_dict(), "git_sha": "flat-sha", "status": "completed"
            }))
            time = np.linspace(0.0, 20.0, 31)
            radius = np.full(len(time), 4.0)
            pd.DataFrame({
                "run_id": f"s{seed}", "time": time, "step": np.arange(len(time)),
                "grain_id": np.ones(len(time), dtype=int), "area": np.pi * radius**2,
                "radius": radius, "perimeter": 2.0 * np.pi * radius,
                "neighbors": np.full(len(time), 6),
            }).to_csv(run_dir / "grain_tracks.csv", index=False)
            pd.DataFrame({
                "time": [0.0, 0.0, 1.0, 1.0],
                "curvature": [0.1] * 4, "normal_velocity": [0.0] * 4,
                "blocked": [1] * 4, "resolved_shear": [0.0] * 4,
                "free_volume_deficit": [0.0] * 4,
                "grain_i": [1] * 4, "grain_j": [2] * 4,
            }).to_csv(run_dir / "boundary_tracks.csv", index=False)
            event_count = 2**temperature_index
            pd.DataFrame({
                "event_type": ["activation_hit"] * event_count,
                "entity_id": ["gb:1"] * event_count,
                "time": np.linspace(0.1, 0.9, event_count),
                "instantaneous_rate": [float(event_count)] * event_count,
                "shear_strain_increment": [0.0] * event_count,
                "volumetric_strain_increment": [0.0] * event_count,
                "barrier_type": ["easy"] * event_count,
                "DeltaG0": [0.2] * event_count,
                "effective_DeltaG": [0.2] * event_count,
                "burgers_vector_b": ["[0.1, 0.0]"] * event_count,
            }).to_csv(run_dir / "events.csv", index=False)
            run_paths.append(str(run_dir))
    campaign = tmp_path / "flat-campaign"
    campaign.mkdir()
    (campaign / "campaign_manifest.json").write_text(json.dumps({
        "runs": run_paths, "status": "completed",
    }))
    summary = analyze_campaign(campaign, bootstrap_samples=5)
    assert summary["n"].isna().all()
    assert (summary["K"] == 0.0).all()
    assert summary["Q_app"].isna().all()
    diagnostics = json.loads((campaign / "mechanism_summary_diagnostics.json").read_text())
    assert all(
        detail["temperature_series_fit"]["kinetically_observable"] is False
        for detail in diagnostics
    )
    assert all(
        detail["temperature_series_fit"]["event_level"]["growth_series_censored"]
        for detail in diagnostics
    )
    assert all(
        detail["temperature_series_fit"]["event_level"]["activation_energy_ev"]
        is not None for detail in diagnostics
    )
