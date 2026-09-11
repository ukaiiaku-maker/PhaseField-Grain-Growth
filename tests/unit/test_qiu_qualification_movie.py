import json

import numpy as np

from render_qiu_qualification_movie import discover, load, select_contact_panels


def test_qiu_movie_discovers_compact_frames_and_loads_required_fields(tmp_path):
    frames = tmp_path / "frames"
    frames.mkdir()
    labels = np.asarray([[0, 1], [1, 0]])
    stress = np.zeros((2, 2, 2, 2)); stress[0, 1] = 2.0
    eigenstrain = np.zeros_like(stress); eigenstrain[0, 1] = -0.5
    path = frames / "frame-0000003.npz"
    np.savez_compressed(
        path, labels=labels, stress=stress, eigenstrain=eigenstrain,
        step=np.asarray(3), time=np.asarray(0.12), grain_count=np.asarray(2),
    )
    assert discover(tmp_path) == [path]
    record = load(path)
    assert record["step"] == 3
    assert record["grain_count"] == 2
    assert np.array_equal(record["stress_xy"], np.full((2, 2), 2.0))
    assert np.array_equal(record["eigenstrain_xy"], np.full((2, 2), -0.5))


def test_qiu_movie_merges_dense_fields_and_prefers_them_at_duplicate_steps(tmp_path):
    frames = tmp_path / "frames"
    fields = tmp_path / "diagnostic_fields"
    frames.mkdir(); fields.mkdir()
    for path in (
        frames / "frame-0000000.npz",
        frames / "frame-0000010.npz",
        fields / "step-0000010-cadence.npz",
        fields / "step-0000011-guard.npz",
    ):
        np.savez_compressed(path, labels=np.zeros((1, 1)))

    assert discover(tmp_path) == [
        frames / "frame-0000000.npz",
        fields / "step-0000010-cadence.npz",
        fields / "step-0000011-guard.npz",
    ]


def test_contact_sheet_does_not_infer_transition_from_sparse_population_loss(tmp_path):
    records = [
        {"step": 0, "grain_count": 800},
        {"step": 200, "grain_count": 700},
        {"step": 400, "grain_count": 610},
        {"step": 500, "grain_count": 580},
    ]
    selected, roles, evidence = select_contact_panels(tmp_path, records)
    assert selected == [0, 2, 3]
    assert roles == ["initial", "midpoint", "terminal"]
    assert evidence == {"selection_basis": "trajectory_without_diagnostic_capture"}


def test_contact_sheet_centers_on_recorded_diagnostic_capture(tmp_path):
    capture = tmp_path / "diagnostic_capture.json"
    capture.write_text(json.dumps({"step": 115, "reasons": ["morphology_guard"]}))
    records = [
        {"step": 0, "grain_count": 800},
        {"step": 100, "grain_count": 790},
        {"step": 110, "grain_count": 780},
        {"step": 120, "grain_count": 700},
    ]
    selected, roles, evidence = select_contact_panels(tmp_path, records)
    assert selected == [1, 2, 3]
    assert roles == ["pre-transition", "transition", "post-transition"]
    assert evidence["selection_basis"] == "diagnostic_capture"
    assert evidence["diagnostic_capture_step"] == 115
    assert evidence["diagnostic_capture_reasons"] == ["morphology_guard"]


def test_explicit_transition_step_overrides_but_preserves_capture_provenance(tmp_path):
    capture = tmp_path / "diagnostic_capture.json"
    capture.write_text(json.dumps({
        "step": 101, "reasons": ["oversensitive_legacy_clipping"],
    }))
    records = [
        {"step": 100}, {"step": 110}, {"step": 120}, {"step": 130},
    ]

    selected, roles, evidence = select_contact_panels(
        tmp_path, records, explicit_transition_step=121,
    )

    assert selected == [1, 2, 3]
    assert roles == ["pre-transition", "transition", "post-transition"]
    assert evidence["selection_basis"] == "explicit_transition_step"
    assert evidence["requested_transition_step"] == 121
    assert evidence["selected_transition_step"] == 120
    assert evidence["preserved_diagnostic_capture_step"] == 101
    assert evidence["preserved_diagnostic_capture_reasons"] == [
        "oversensitive_legacy_clipping"
    ]
