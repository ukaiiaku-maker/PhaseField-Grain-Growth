from __future__ import annotations

import numpy as np

from grain_growth_pf.analysis.kinetic_resistance import fit_inverse_growth_rate


def test_inverse_growth_rate_fit_recovers_labelled_law():
    radius = np.linspace(2.0, 8.0, 25)
    expected_a = 1.75
    expected_b = 0.4
    growth_rate = 1.0 / (expected_a * radius + expected_b)

    fit = fit_inverse_growth_rate(radius, growth_rate)

    assert np.isclose(fit.a, expected_a)
    assert np.isclose(fit.b, expected_b)
    assert np.isclose(fit.r_squared, 1.0)
