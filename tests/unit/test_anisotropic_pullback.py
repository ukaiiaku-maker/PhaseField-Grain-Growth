import numpy as np
import pytest

from grain_growth_pf.mechanics.anisotropy import BoundaryLaw, LADDER
from grain_growth_pf.mechanics.network_pullback import triangle_energy_pullback


@pytest.mark.parametrize('values', [
    [[.8, .2, .3], [.2, .8, .7]],
    [[.9, .1, .05], [.07, .83, .13], [.03, .07, .82]],
])
def test_phase_pullback_matches_independent_energy_perturbations(values):
    values = np.asarray(values)
    coordinates = np.array([[.2, -.1], [1.4, .2], [.1, 1.1]])
    orientations = np.array([.1, .3, .7])[:len(values)]
    law = BoundaryLaw()
    energy, gradient = triangle_energy_pullback(values, coordinates, orientations, law)
    assert energy > 0
    np.testing.assert_allclose(gradient.sum(axis=0), 0, atol=1e-13)
    for eps in (1e-4, 1e-5, 1e-6):
        measured = np.zeros_like(values)
        for phase, vertex in np.ndindex(values.shape):
            plus, minus = values.copy(), values.copy()
            plus[phase, vertex] += eps
            minus[phase, vertex] -= eps
            measured[phase, vertex] = (
                triangle_energy_pullback(plus, coordinates, orientations, law)[0]
                - triangle_energy_pullback(minus, coordinates, orientations, law)[0])/(2*eps)
        assert np.linalg.norm(measured-gradient)/np.linalg.norm(gradient) < 1e-4


def test_isotropic_reference_subtraction_is_exact():
    values = np.array([[.8, .2, .3], [.2, .8, .7]])
    coordinates = [[0, 0], [1, 0], [0, 1]]
    reference = BoundaryLaw(LADDER['A0_ISOTROPIC'])
    _, original = triangle_energy_pullback(values, coordinates, [.1, .3], reference)
    _, nested = triangle_energy_pullback(values, coordinates, [.7, 1.2], reference)
    np.testing.assert_array_equal(original-nested, 0)
