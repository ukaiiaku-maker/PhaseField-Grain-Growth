import numpy as np

from grain_growth_pf.encounters.gb_area import point_defect_requirement
from grain_growth_pf.encounters.geometric_hazard import GeometricEncounterClock
from grain_growth_pf.encounters.swept_volume import swept_measure
from grain_growth_pf.encounters.tj_path import periodic_path_increment
from grain_growth_pf.entities.tracker import EntityTracker


def test_controlled_boundary_length_and_no_pixel_churn_hazard():
    labels = np.zeros((20, 30), dtype=int)
    labels[:, 15:] = 1
    tracker = EntityTracker(np.array([0.0, 0.1]), dx=0.5, domain_length=100, periodic=False)
    first = tracker.update(labels)
    segment = next(iter(first.boundaries.values()))
    assert np.isclose(segment.length, 20 * 0.5)
    segment.encounter_state["hazard"] = 1.2
    second = tracker.update(labels.copy())
    assert next(iter(second.boundaries.values())).encounter_state["hazard"] == 1.2


def test_path_swept_area_and_quota_exact():
    assert np.isclose(periodic_path_increment((9, 5), (1, 5), (10, 10)), 2)
    assert swept_measure(8, 0.5, physics_dimension=2) == 4
    assert swept_measure(8, 0.5, physics_dimension=3, out_of_plane_thickness=2) == 8
    assert point_defect_requirement(-4, 0.2, 0.05) == 16


def test_geometric_encounter_stops_at_physical_state_change():
    stopped = GeometricEncounterClock(2.0, np.random.default_rng(201))
    stopped.threshold = 0.5
    events = stopped.advance(10.0, maximum_events=1)
    assert len(events) == 1
    assert events[0].overshoot == 0.0
    assert np.isclose(events[0].measure, 0.25)
    assert np.isclose(stopped.cumulative_hazard, 0.5)
    assert np.isclose(stopped.total_measure, 0.25)

    continuous = GeometricEncounterClock(2.0, np.random.default_rng(201))
    continuous.threshold = 0.5
    events = continuous.advance(10.0)
    assert len(events) > 1
    assert np.isclose(continuous.cumulative_hazard, 20.0)
    assert np.isclose(continuous.total_measure, 10.0)
