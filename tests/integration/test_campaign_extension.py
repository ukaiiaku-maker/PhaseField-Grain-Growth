import json
from pathlib import Path

import numpy as np

from grain_growth_pf.campaign import extend_campaign
from grain_growth_pf.config import ModelConfig, PFConfig
from grain_growth_pf.simulation import EventResolvedSimulation


def test_campaign_extension_is_exact_and_preserves_source(tmp_path):
    common = dict(
        regime="B0",
        seed=81,
        pf=PFConfig(
            shape=(20, 20), interface_width=3, time_step=0.01,
            intrinsic_mobility=0.2, adaptive_stepping=True,
        ),
        output_cadence=1,
        termination_grains=1,
        parameters={"initial_grains": 7},
    )
    source_config = ModelConfig(**common, max_steps=2)
    source_run = tmp_path / "source-run"
    EventResolvedSimulation(source_config, source_run).run()
    source_eta = np.load(source_run / "checkpoint.npz")["eta"].copy()
    source_tracks = (source_run / "grain_tracks.csv").read_bytes()
    source_campaign = tmp_path / "source-campaign"
    source_campaign.mkdir()
    (source_campaign / "campaign_manifest.json").write_text(json.dumps({
        "runs": [str(source_run)], "status": "completed"
    }))

    continuous_config = ModelConfig(**common, max_steps=4)
    continuous = EventResolvedSimulation(continuous_config, tmp_path / "continuous")
    continuous.run()
    extended_campaign = extend_campaign(
        [source_campaign], max_steps=4, termination_grains=1,
        root=tmp_path / "extensions", processes=1,
    )
    campaign_manifest = json.loads((extended_campaign / "campaign_manifest.json").read_text())
    extended_run = campaign_manifest["runs"][0]
    extended_eta = np.load(Path(extended_run) / "checkpoint.npz")["eta"]

    assert np.array_equal(extended_eta, continuous.solver.eta)
    assert np.array_equal(np.load(source_run / "checkpoint.npz")["eta"], source_eta)
    assert (source_run / "grain_tracks.csv").read_bytes() == source_tracks
    assert campaign_manifest["status"] == "completed"
    run_manifest = json.loads((Path(extended_run) / "manifest.json").read_text())
    assert run_manifest["restart_provenance"]["checkpoint_step"] == 2
    assert run_manifest["steps_completed"] == 4
