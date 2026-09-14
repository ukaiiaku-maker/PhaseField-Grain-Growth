from __future__ import annotations

import json

import numpy as np

from grain_growth_pf.config import PFConfig
from grain_growth_pf.pf.solver import MultiphaseFieldSolver
from scripts.recover_anisotropic_postfix import _load_state, _save_state, SOURCE_COMMIT


def test_recovery_checkpoint_round_trip_contains_restart_state(tmp_path):
    eta = np.zeros((2, 4, 4), dtype=float)
    eta[0, :, :2] = 1.0
    eta[1, :, 2:] = 1.0
    cfg = PFConfig(shape=(4, 4), time_step=1e-3)
    solver = MultiphaseFieldSolver(eta, cfg, orientations=np.array([0.0, 0.2]))
    solver.step(compute_energy=False)
    metadata = {
        "schema": "anisotropic-postfix-checkpoint-v1",
        "source_commit": SOURCE_COMMIT,
        "case": "unit",
        "stage": "continuous",
        "accepted_step": solver.step_number,
        "physical_time": solver.time,
        "identities": {"input": "test"},
    }
    path = tmp_path / "checkpoint.npz"
    _save_state(path, solver, metadata)

    restored = MultiphaseFieldSolver(eta, cfg, orientations=np.array([0.0, 0.2]))
    loaded = _load_state(path, restored)

    assert {key: loaded[key] for key in metadata} == metadata
    assert loaded["solver_config"]["shape"] == [4, 4]
    assert np.array_equal(restored.eta, solver.eta)
    assert np.array_equal(restored.active_phases, solver.active_phases)
    assert np.array_equal(restored.mobility_scale, solver.mobility_scale)
    assert restored.step_number == solver.step_number
    assert restored.time == solver.time


def test_atomic_checkpoint_metadata_is_embedded(tmp_path):
    eta = np.stack((np.ones((3, 3)), np.zeros((3, 3))))
    solver = MultiphaseFieldSolver(
        eta, PFConfig(shape=(3, 3)), orientations=np.array([0.0, 0.2])
    )
    metadata = {
        "schema": "anisotropic-postfix-checkpoint-v1",
        "source_commit": SOURCE_COMMIT,
        "case": "unit",
        "stage": "midpoint",
        "accepted_step": 0,
        "physical_time": 0.0,
        "identities": {},
    }
    path = tmp_path / "state.npz"
    _save_state(path, solver, metadata)
    with np.load(path, allow_pickle=False) as archive:
        embedded = json.loads(str(archive["checkpoint_state_json"]))
        assert {key: embedded[key] for key in metadata} == metadata
        assert embedded["solver_config"]["shape"] == [3, 3]
    assert not list(tmp_path.glob("*.tmp"))
