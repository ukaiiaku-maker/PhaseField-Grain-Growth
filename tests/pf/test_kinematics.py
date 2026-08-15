import numpy as np

from grain_growth_pf.pf.geometry import circular_grain, planar_interface
from grain_growth_pf.pf.kinematics import interface_kinematics


def _interface_points(phase):
    return np.argwhere((phase > 0.35) & (phase < 0.65)).astype(float)


def test_diffuse_circle_kinematics_has_capillary_sign_and_magnitude():
    previous = circular_grain((96, 96), 20.1, 4.0)[1]
    current = circular_grain((96, 96), 20.0, 4.0)[1]
    curvature, velocity, normal = interface_kinematics(
        current, previous, _interface_points(current), 0.1, 1.0
    )
    assert np.isclose(curvature, -1 / 20, rtol=0.2)
    assert np.isclose(velocity, -1.0, rtol=0.12)
    assert np.linalg.norm(normal) < 0.2


def test_translated_planar_interface_has_zero_curvature():
    previous = planar_interface((64, 64), 4.0, angle=0)[1]
    current = np.roll(previous, 1, axis=1)
    curvature, velocity, _ = interface_kinematics(
        current, previous, _interface_points(current), 1.0, 1.0, periodic=True
    )
    assert abs(curvature) < 1e-12
    assert np.isfinite(velocity)
