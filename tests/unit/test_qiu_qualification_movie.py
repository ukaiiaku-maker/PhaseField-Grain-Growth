import numpy as np

from render_qiu_qualification_movie import discover, load


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
