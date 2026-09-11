from pathlib import Path

import pytest

from grain_growth_pf.config import ModelConfig, PFConfig
from continue_qiu_legacy_forensics import continuation_config
from run_qiu_fft_v2_qualification import apply_run_overrides


def _config() -> ModelConfig:
    return ModelConfig(
        regime="FFT_EIGENSTRAIN_V2", seed=5101,
        pf=PFConfig(shape=(16, 16), time_step=0.04),
        mechanics_backend="fft_eigenstrain_v2",
        parameters={"external_delta_eta_target": 0.02},
    )


def test_refinement_overrides_actual_timestep_and_external_target(tmp_path):
    updated = apply_run_overrides(
        _config(), tmp_path / "initial.npz", time_step=0.02, target=0.01,
    )

    assert updated.pf.time_step == 0.02
    assert updated.parameters["external_delta_eta_target"] == 0.01
    assert updated.parameters["initial_state_file"] == str((tmp_path / "initial.npz").resolve())


@pytest.mark.parametrize("name,value", [("time_step", 0.0), ("target", -0.01)])
def test_nonpositive_refinement_override_is_rejected(name: str, value: float):
    with pytest.raises(ValueError):
        apply_run_overrides(_config(), Path("initial.npz"), **{name: value})


def test_legacy_continuation_changes_diagnostics_but_not_backend_or_timestep(tmp_path):
    legacy = ModelConfig(
        regime="QIU", seed=5101,
        pf=PFConfig(shape=(16, 16), time_step=0.04),
        mechanics_backend="qiu_full_field", output_cadence=200,
        max_steps=100000, termination_grains=100,
        parameters={"checkpoint_cadence": 1000, "easy_beta": 0.35},
    )
    initial = tmp_path / "initial.npz"
    continued = continuation_config(
        {"config": legacy.to_dict()}, initial, end_step=10246,
    )

    assert continued.mechanics_backend == legacy.mechanics_backend
    assert continued.pf == legacy.pf
    assert continued.parameters["easy_beta"] == legacy.parameters["easy_beta"]
    assert continued.parameters["qiu_guard_terminate"] is False
    assert continued.parameters["qiu_guard_clip_warmup_steps"] == 20
    assert continued.parameters["qiu_diagnostic_field_cadence"] == 10
    assert continued.parameters["checkpoint_cadence"] == 500
    assert continued.max_steps == 10246
