import numpy as np

from audit_qiu_transition_fragment import domain_measurements, relative_l2


def test_domain_measurements_exposes_mismatched_gradient_singularity():
    phase = np.zeros((7, 7), dtype=float)
    previous = np.zeros_like(phase)
    partner = np.zeros_like(phase)
    # Radius-one opposing samples cancel exactly, while radius two is nonzero.
    phase[3, 3] = 0.5
    phase[1, 3] = 0.0
    phase[5, 3] = 1.0
    partner[3, 3] = 0.5
    result = domain_measurements(
        phase, previous, partner, np.asarray([[3.0, 3.0]]),
        0.04, 1.0, periodic=True,
    )
    point = result["points"][0]
    assert point["valid"] is True
    assert point["gradient_radius1"] == 0.0
    assert point["gradient_radius2"] == 0.25
    assert point["point_velocity"] == 1.25e15


def test_relative_l2_uses_reference_scale():
    reference = np.asarray([3.0, 4.0])
    assert relative_l2(np.asarray([0.3, 0.4]), reference) == 0.1
