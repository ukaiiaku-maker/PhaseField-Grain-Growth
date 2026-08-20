import numpy as np

from grain_growth_pf.entities.arclength_tracker import ArclengthEntityTracker


def _inclusion_labels(shift_y=0, shift_x=0):
    labels = np.zeros((32, 32), dtype=int)
    y0, y1 = 8 + shift_y, 24 + shift_y
    x0, x1 = 8 + shift_x, 24 + shift_x
    labels[y0:y1, x0:x1] = 1
    return labels


def test_arclength_domains_are_connected_and_physically_bounded():
    tracker = ArclengthEntityTracker(
        np.array([0.0, 0.3]), dx=1.0, domain_length=6.0, periodic=False
    )
    snapshot = tracker.update(_inclusion_labels())
    segments = [
        segment for segment in snapshot.boundaries.values()
        if {segment.grain_i, segment.grain_j} == {0, 1}
    ]
    assert len(segments) > 4
    for segment in segments:
        if len(segment.points) > 1:
            increments = np.abs(np.diff(segment.points, axis=0))
            assert np.all(np.max(increments, axis=1) <= 1.0)
        assert segment.length <= 7.0


def test_arclength_domain_state_persists_under_small_boundary_motion():
    tracker = ArclengthEntityTracker(
        np.array([0.0, 0.3]), dx=1.0, domain_length=8.0, periodic=False
    )
    first = tracker.update(_inclusion_labels())
    key = sorted(first.boundaries)[0]
    first.boundaries[key].shear_incompatibility = 4.25
    second = tracker.update(_inclusion_labels(shift_x=1))
    assert key in second.boundaries
    assert second.boundaries[key].shear_incompatibility == 4.25
