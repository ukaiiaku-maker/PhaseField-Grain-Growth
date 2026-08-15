import csv

import numpy as np

from grain_growth_pf.analysis.activation_energy import fit_activation_energy
from grain_growth_pf.analysis.analytical_models import asymptotic_exponent, intrinsic_radius, poisson_activity, series_activity
from grain_growth_pf.analysis.growth_law import fit_growth_law
from grain_growth_pf.disconnections.mode import K_B_EV
from grain_growth_pf.io.event_ledger import EVENT_FIELDS, EventLedger


def test_growth_and_activation_recover_inputs():
    time = np.linspace(0, 100, 101)
    radius = (4**3 + 0.7 * time) ** (1 / 3)
    fit = fit_growth_law(time, radius)
    assert abs(fit.exponent - 3) < 1e-3
    assert fit.r_squared > 0.999999
    # A small relative radius change must not bias the fit toward n=1.
    narrow_time = np.linspace(0, 20, 80)
    narrow_radius = (100 + 0.3 * narrow_time) ** 0.5
    assert abs(fit_growth_law(narrow_time, narrow_radius).exponent - 2) < 1e-3
    temps = np.array([700, 800, 900, 1050, 1200])
    q = 0.65
    coefficients = 3e5 * np.exp(-q / (K_B_EV * temps))
    activation = fit_activation_energy(temps, coefficients)
    assert abs(activation.activation_energy_ev - q) < 1e-10


def test_analytical_limits():
    t = np.arange(4.0)
    assert np.allclose(intrinsic_radius(t, 2, 0.5) ** 2, 4 + t)
    assert np.isclose(poisson_activity(1, 2), 1 - np.exp(-2))
    assert np.isclose(series_activity(0.5, 0.25), 1 / 6)
    assert asymptotic_exponent(1, 3) == 5


def test_event_ledger_schema(tmp_path):
    target = tmp_path / "events.csv"
    with EventLedger(target) as ledger:
        ledger.write({"run_id": "x", "event_id": 1, "normal_step_h": 0.2})
    with target.open() as handle:
        row = next(csv.DictReader(handle))
    assert tuple(row) == EVENT_FIELDS
    assert row["normal_step_h"] == "0.2"
