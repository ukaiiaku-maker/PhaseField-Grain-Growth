import json

import pytest

from grain_growth_pf.campaign import compose_completed_campaigns, enumerate_campaign


def test_campaign_regime_overrides_merge_parameters_and_pf():
    spec = {
        "base": {
            "pf": {"shape": [24, 24], "temperature": 800},
            "parameters": {"initial_grains": 8, "shared": 1},
            "max_steps": 10,
        },
        "regimes": {
            "A": {},
            "B": {"parameters": {"shared": 2, "barrier_ev": 0.5},
                  "pf": {"intrinsic_mobility": 3}, "max_steps": 20},
        },
        "temperatures": [700, 900],
        "seeds": [4],
    }
    configs = enumerate_campaign(spec)
    assert len(configs) == 4
    b = next(config for config in configs if config.regime == "B" and config.pf.temperature == 900)
    assert b.parameters == {"initial_grains": 8, "shared": 2, "barrier_ev": 0.5}
    assert b.pf.intrinsic_mobility == 3
    assert b.max_steps == 20


def test_composite_campaign_accepts_only_unique_completed_runs(tmp_path):
    sources = []
    for source_index in range(2):
        source = tmp_path / f"source-{source_index}"
        source.mkdir()
        runs = []
        for seed, status in ((source_index + 1, "completed"), (20 + source_index, "running")):
            run = source / f"run-{seed}"
            run.mkdir()
            (run / "manifest.json").write_text(json.dumps({
                "status": status,
                "config": {
                    "regime": "R", "seed": seed,
                    "pf": {"temperature": 900.0},
                },
            }))
            runs.append(str(run))
        (source / "campaign_manifest.json").write_text(json.dumps({"runs": runs}))
        sources.append(source)

    composite = compose_completed_campaigns(sources, tmp_path / "composites", expected_runs=2)
    manifest = json.loads((composite / "campaign_manifest.json").read_text())
    assert manifest["status"] == "completed"
    assert manifest["runs_total"] == 2
    assert len(manifest["runs"]) == 2
    with pytest.raises(ValueError, match="expected 3"):
        compose_completed_campaigns(sources, tmp_path / "other", expected_runs=3)
