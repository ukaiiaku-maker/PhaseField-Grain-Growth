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


def test_boundary_domain_split_merge_and_tj_reconnection_retire_state():
    tracker = EntityTracker(
        np.array([0.0, 0.2, 0.4]), domain_length=100, periodic=False
    )
    initial = tracker.update(three_grain_labels())
    pair_key = next(key for key in initial.boundaries if "0-1" in key)
    initial.boundaries[pair_key].free_volume_deficit = 4.0
    old_tj = next(iter(initial.triple_junctions.values()))
    old_tj.add_burgers(np.array([0.3, -0.2]))

    tracker.domain_length = 2.0
    split = tracker.update(three_grain_labels())
    split_keys = sorted(key for key in split.boundaries if "0-1" in key)
    assert len(split_keys) > 1
    assert split.boundaries[split_keys[0]].free_volume_deficit == 4.0
    assert all(
        split.boundaries[key].free_volume_deficit == 0.0 for key in split_keys[1:]
    )
    for tj in split.triple_junctions.values():
        assert len(tj.adjoining_boundaries) == 3
        adjoining_pairs = {
            tuple(sorted((split.boundaries[key].grain_i, split.boundaries[key].grain_j)))
            for key in tj.adjoining_boundaries
        }
        assert adjoining_pairs == {(0, 1), (0, 2), (1, 2)}

    tracker.domain_length = 100.0
    merged = tracker.update(three_grain_labels())
    merged_keys = [key for key in merged.boundaries if "0-1" in key]
    assert merged_keys == [pair_key]
    assert merged.boundaries[pair_key].free_volume_deficit == 4.0

    without_third = three_grain_labels()
    without_third[without_third == 2] = 0
    assert not tracker.update(without_third).triple_junctions
    reconnected = tracker.update(np.roll(three_grain_labels(), 1, axis=0))
    new_tj = next(iter(reconnected.triple_junctions.values()))
    assert np.array_equal(new_tj.residual_burgers, [0.0, 0.0])
    assert new_tj.travel_distance == 0.0
