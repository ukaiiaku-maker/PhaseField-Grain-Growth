from pathlib import Path

import pytest

from grain_growth_pf.config import ModelConfig, PFConfig
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
