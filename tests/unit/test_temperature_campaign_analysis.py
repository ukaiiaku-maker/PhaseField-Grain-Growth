import json

import numpy as np
import pandas as pd

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
    (campaign / "campaign_manifest.json").write_text(json.dumps({"runs": run_paths}))

    summary = analyze_campaign(campaign, bootstrap_samples=5)
    assert np.allclose(summary["n"], 2.0, atol=1e-6)
    assert np.allclose(summary.sort_values("temperature")["K"], coefficients, rtol=1e-6)
    assert np.allclose(summary["Q_app"], activation_energy, atol=1e-6)
    diagnostics = json.loads((campaign / "mechanism_summary_diagnostics.json").read_text())
    assert all("temperature_series_fit" in detail for detail in diagnostics)
