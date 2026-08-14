import numpy as np

from grain_growth_pf.entities.tracker import EntityTracker


def three_grain_labels():
    labels = np.zeros((12, 12), dtype=int)
    labels[:, 6:] = 1
    labels[6:, 3:9] = 2
    return labels


def test_entity_birth_motion_and_no_state_transfer():
    tracker = EntityTracker(np.array([0.0, 0.2, 0.4, 0.6]), domain_length=100, periodic=False)
    first = tracker.update(three_grain_labels())
    key = next(k for k in first.boundaries if "0-1" in k)
    first.boundaries[key].shear_incompatibility = 7.0
    moved = np.roll(three_grain_labels(), 1, axis=1)
    second = tracker.update(moved)
    assert second.boundaries[key].shear_incompatibility == 7.0
    assert any(set(t.grain_ids) == {0, 1, 2} for t in second.triple_junctions.values())
    # Remove grain 1, then introduce unrelated grain 3: retired state must not transfer.
    gone = moved.copy(); gone[gone == 1] = 0
    tracker.update(gone)
    born = gone.copy(); born[:, 8:] = 3
    final = tracker.update(born)
    assert not any("0-1" in k for k in final.boundaries)
    assert all(b.shear_incompatibility == 0 for b in final.boundaries.values())


def test_tj_path_and_burgers_persist():
    tracker = EntityTracker(np.array([0.0, 0.2, 0.4]), domain_length=100, periodic=False)
    snap = tracker.update(three_grain_labels())
    tj = next(iter(snap.triple_junctions.values()))
    tj.add_burgers(np.array([1.0, -0.5]))
    moved = np.roll(three_grain_labels(), 1, axis=0)
    snap2 = tracker.update(moved)
    tj2 = snap2.triple_junctions[tj.entity_id]
    assert tj2.travel_distance > 0
    assert np.array_equal(tj2.residual_burgers, [1.0, -0.5])

