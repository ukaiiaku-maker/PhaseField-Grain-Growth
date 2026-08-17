import json

import pytest

from grain_growth_pf.campaign import resume_campaign
from grain_growth_pf.config import ModelConfig, PFConfig
from grain_growth_pf.simulation import EventResolvedSimulation


def test_resume_campaign_continues_checkpoints_and_starts_queued_runs(tmp_path):
    config = ModelConfig(
        regime="B0", seed=81,
        pf=PFConfig(
            shape=(18, 18), interface_width=3, time_step=0.01,
            intrinsic_mobility=0.2, adaptive_stepping=True,
        ),
        output_cadence=1, max_steps=3, termination_grains=1,
        parameters={"initial_grains": 6},
    )
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    resumed_run = campaign / "B0-T900-s81-existing"
    queued_run = campaign / "B0-T900-s82-queued"
    simulation = EventResolvedSimulation(config, resumed_run, code_sha="source-sha")
    simulation.solver.step()
    simulation.snapshot = simulation.tracker.update(simulation.solver.labels)
    simulation._update_physics()
    simulation._write_tracks()
    simulation._save_checkpoint()
    simulation.ledger.close()
    simulation.track_handle.close()
    simulation.boundary_handle.close()

    specification = {
        "temperatures": [900.0],
        "seeds": [81, 82],
        "base": config.to_dict(),
        "regimes": {"B0": {}},
    }
    (campaign / "campaign_manifest.json").write_text(json.dumps({
        "source_spec": str(tmp_path / "spec.yaml"),
        "specification": specification,
        "initial_condition_files": {},
        "runs": [str(resumed_run), str(queued_run)],
        "status": "running",
    }))

    resume_campaign(campaign, processes=1)

    manifest = json.loads((campaign / "campaign_manifest.json").read_text())
    assert manifest["status"] == "completed"
    assert manifest["restart_history"][-1]["resumed_runs"] == [str(resumed_run)]
    assert manifest["restart_history"][-1]["started_queued_runs"] == [str(queued_run)]
    for run in (resumed_run, queued_run):
        run_manifest = json.loads((run / "manifest.json").read_text())
        assert run_manifest["status"] == "completed"
        assert run_manifest["steps_completed"] == 3
    resumed_manifest = json.loads((resumed_run / "manifest.json").read_text())
    assert resumed_manifest["restart_history"][-1]["checkpoint_step"] == 1
    assert resumed_manifest["restart_history"][-1]["source_git_sha"] == "source-sha"


def test_keyboard_interrupt_is_not_reported_as_completion(tmp_path):
    config = ModelConfig(
        regime="interrupt", seed=90,
        pf=PFConfig(shape=(14, 14), interface_width=3, time_step=0.01),
        output_cadence=1, max_steps=2, termination_grains=1,
        parameters={"initial_grains": 4},
    )
    run = tmp_path / "interrupted"
    simulation = EventResolvedSimulation(config, run)

    def interrupt():
        raise KeyboardInterrupt

    simulation.solver.step = interrupt
    with pytest.raises(KeyboardInterrupt):
        simulation.run()
    manifest = json.loads((run / "manifest.json").read_text())
    assert manifest["status"] == "failed"
    assert manifest["failure"].startswith("KeyboardInterrupt")
