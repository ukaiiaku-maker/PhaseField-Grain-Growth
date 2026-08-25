from __future__ import annotations

import numpy as np

from grain_growth_pf.analysis.long_time import (
    differential_effective_exponent,
    grain_size_moments,
    profile_growth_window,
)


def test_grain_size_moments_match_definitions():
    areas = np.pi * np.asarray([1.0, 4.0, 9.0])
    moments = grain_size_moments(areas)
    assert np.isclose(moments.number_mean, 4.0)
    assert np.isclose(moments.population, 2.0 * np.sqrt(14.0 / 3.0))
    assert np.isclose(moments.area_weighted, 72.0 / 14.0)


def test_profile_growth_window_recovers_known_cubic_law():
    time = np.linspace(0.0, 100.0, 101)
    size = (8.0**3 + 2.5 * time) ** (1.0 / 3.0)
    fit = profile_growth_window(time, size)
    assert np.isclose(fit.n_best, 3.0, atol=0.1)
    assert np.isclose(fit.coefficient, 2.5, rtol=2e-3)
    assert fit.normalized_rmse < 1e-10
    assert not fit.hit_search_bound


def test_profile_growth_window_extends_exponent_search_and_flags_bound():
    time = np.linspace(0.0, 1.0, 101)
    size = (1.0 + time) ** (1.0 / 50.0)
    fit = profile_growth_window(time, size, n_max=20.0, extended_n_max=50.0)
    assert fit.n_best >= 49.0
    assert fit.hit_search_bound


def test_differential_effective_exponent_is_cubic_away_from_edges():
    time = np.linspace(0.0, 100.0, 201)
    size = (8.0**3 + 2.5 * time) ** (1.0 / 3.0)
    rate, exponent = differential_effective_exponent(time, size, window_length=21)
    assert np.all(rate[20:-20] > 0.0)
    np.testing.assert_allclose(np.nanmedian(exponent[20:-20]), 3.0, atol=0.08)
