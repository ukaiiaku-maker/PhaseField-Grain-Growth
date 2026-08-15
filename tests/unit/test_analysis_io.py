import csv
import json

import numpy as np
import pandas as pd

from grain_growth_pf.analysis.activation_energy import fit_activation_energy
from grain_growth_pf.analysis.analytical_models import asymptotic_exponent, intrinsic_radius, poisson_activity, series_activity
from grain_growth_pf.analysis.growth_law import (
    fit_common_exponent,
    fit_growth_law,
    fit_growth_law_fixed_exponent,
    scan_growth_exponent,
)
from grain_growth_pf.analysis.grain_tracks import ensemble_radius
from grain_growth_pf.disconnections.mode import K_B_EV
from grain_growth_pf.io.event_ledger import EVENT_FIELDS, EventLedger
from grain_growth_pf.io.provenance import write_manifest


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

    profile = scan_growth_exponent(time, radius, np.linspace(1.0, 4.0, 301))
    assert abs(profile.exponents[np.argmin(profile.normalized_rmse)] - 3) < 0.02
    fixed = fit_growth_law_fixed_exponent(time, radius, 3.0)
    assert abs(fixed.coefficient - 0.7) < 1e-10
    common = fit_common_exponent(
        [time, time, time],
        [(4**3 + coefficient * time) ** (1 / 3) for coefficient in (0.2, 0.7, 1.4)],
    )
    assert abs(common.exponent - 3.0) < 1e-7
    assert np.allclose(common.coefficients, [0.2, 0.7, 1.4], rtol=1e-8)
    temps = np.array([700, 800, 900, 1050, 1200])
    q = 0.65
    coefficients = 3e5 * np.exp(-q / (K_B_EV * temps))
    activation = fit_activation_energy(temps, coefficients)
    assert abs(activation.activation_energy_ev - q) < 1e-10


def test_ensemble_radius_reports_independent_size_measures():
    tracks = pd.DataFrame({
        "run_id": ["x", "x"], "time": [0.0, 0.0], "step": [0, 0],
        "grain_id": [1, 2], "area": [np.pi, 9 * np.pi],
        "radius": [1.0, 3.0], "perimeter": [2 * np.pi, 8 * np.pi],
    })
    row = ensemble_radius(tracks).iloc[0]
    assert np.isclose(row["R_A"], np.sqrt(5.0))
    assert np.isclose(row["R_mean"], 2.0)
    assert np.isclose(row["R_median"], 2.0)
    assert np.isclose(row["R_rms"], np.sqrt(5.0))
    assert np.isclose(row["R_perimeter"], 2.5)


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


def test_manifest_can_pin_launch_revision(tmp_path):
    target = tmp_path / "manifest.json"
    write_manifest(target, {"regime": "test"}, "running", code_sha="launch-sha")
    write_manifest(target, {"regime": "test"}, "completed", code_sha="launch-sha")
    manifest = json.loads(target.read_text())
    assert manifest["git_sha"] == "launch-sha"
    assert manifest["status"] == "completed"
