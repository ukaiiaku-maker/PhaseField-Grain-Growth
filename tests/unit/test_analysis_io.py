import csv
import json

import numpy as np
import pandas as pd

from grain_growth_pf.analysis.activation_energy import (
    fit_activation_energy,
    local_activation_energies,
)
from grain_growth_pf.analysis.analytical_models import (
    asymptotic_exponent,
    crossover_radius_prediction,
    fit_crossover_growth,
    intrinsic_radius,
    poisson_activity,
    series_activity,
)
from grain_growth_pf.analysis.growth_law import (
    fit_common_exponent,
    fit_growth_law,
    fit_growth_law_fixed_exponent,
    scan_growth_exponent,
)
from grain_growth_pf.analysis.grain_tracks import ensemble_radius
from grain_growth_pf.analysis.campaign import (
    _boundary_metrics,
    _burst_size_ccdf,
    _event_diagnostics,
    _event_rate_observation,
    _fit_window,
    _spatial_motion_correlation,
    _trajectory_distributions,
)
from grain_growth_pf.disconnections.mode import K_B_EV
from grain_growth_pf.io.event_ledger import (
    EVENT_FIELDS,
    EventLedger,
    event_ledger_has_rows,
    event_ledger_path,
    read_event_ledger,
)
from grain_growth_pf.io.provenance import write_manifest
from grain_growth_pf.analysis.jerkiness import jerkiness_metrics
from grain_growth_pf.analysis.plots import (
    _arrhenius_figure,
    _local_exponent,
    _tj_failure_figure,
)


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
    _, local_q = local_activation_energies(temps, coefficients)
    assert np.allclose(local_q, q)

    local = _local_exponent(time, radius, half_window=20)
    assert np.allclose(local[np.isfinite(local)], 3.0, atol=0.05)


def test_arrhenius_plot_includes_global_and_local_diagnostics(tmp_path):
    temperatures = np.asarray([800.0, 900.0, 1000.0, 1100.0])
    coefficients = np.asarray([0.5, 0.8, 1.1, 1.4])
    summary = pd.DataFrame({
        "temperature": temperatures, "K": coefficients,
        "K_ci": 0.05 * coefficients,
    })
    detail = {
        "temperatures": temperatures.tolist(),
        "activation_energy_ev": 0.42,
        "activation_energy_95pct": [0.38, 0.46],
        "local_activation_midpoint_temperature": [850.0, 950.0, 1050.0],
        "local_activation_energy_ev": [0.40, 0.42, 0.44],
        "event_level": {
            "rates": [1.0, 2.0, 4.0, 8.0],
            "activation_energy_ev": 0.55,
            "activation_energy_95pct": [0.50, 0.60],
            "local_activation_midpoint_temperature": [850.0, 950.0, 1050.0],
            "local_activation_energy_ev": [0.50, 0.55, 0.60],
        },
    }
    _arrhenius_figure("synthetic", summary, detail, tmp_path / "arrhenius")
    assert (tmp_path / "arrhenius.png").exists()
    assert (tmp_path / "arrhenius.pdf").exists()


def test_tj_failure_plot_separates_bare_and_residual_adjusted_barriers(tmp_path):
    pd.DataFrame({
        "event_type": ["tj_compatibility_failure"] * 3 + ["disconnection_mode"],
        "entity_id": ["tj:1", "tj:1", "tj:2", "gb:1"],
        "barrier_type": ["easy", "easy", "rare", "easy"],
        "DeltaG0": [0.2, 0.2, 0.8, 0.2],
        "effective_DeltaG": [0.25, 0.18, 1.0, 0.2],
        "burgers_vector_b": ["[0.1, 0.0]", "[-0.1, 0.0]", "[0.0, 0.2]", ""],
    }).to_csv(tmp_path / "events.csv", index=False)
    target = tmp_path / "tj-failures"
    assert _tj_failure_figure([tmp_path], target)
    assert target.with_suffix(".png").exists()
    assert target.with_suffix(".pdf").exists()


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

    resumed = pd.concat([
        tracks.assign(run_id="restart", time=1.0, step=1),
        tracks.assign(run_id="original"),
    ], ignore_index=True)
    resumed_radius = ensemble_radius(resumed)
    assert resumed_radius["time"].tolist() == [0.0, 1.0]
    assert resumed_radius["grain_count"].tolist() == [2, 2]


def test_topology_window_uses_broad_post_equilibration_interval():
    shallow = np.linspace(200, 110, 91)
    start, end, reason = _fit_window(shallow)
    assert shallow[start] <= 190
    assert end == len(shallow)
    assert reason == "five_pct_loss_to_available_end"


def test_jerkiness_reports_motion_concentration_and_stationarity():
    metrics = jerkiness_metrics(
        np.arange(6.0), np.array([0.0, 0.0, 0.0, 10.0, 10.0, 10.0]),
        events=np.array([0, 0, 3, 0, 0]),
    )
    assert np.isclose(metrics["stationary_fraction"], 0.8)
    assert np.isclose(metrics["motion_top_1pct"], 1.0)
    assert np.isclose(metrics["motion_top_5pct"], 1.0)
    assert np.isclose(metrics["motion_top_10pct"], 1.0)
    assert metrics["Fano"] > 1.0

    deep = np.linspace(200, 60, 141)
    start, end, reason = _fit_window(deep)
    assert deep[start] <= 190
    assert deep[start - 1] > 190
    assert end == len(deep)
    assert reason == "five_pct_loss_to_available_end"


def test_boundary_reverse_metric_filters_inactive_diffuse_jitter(tmp_path):
    values = pd.DataFrame({
        "curvature": [2.0, 1.8, 0.01, -0.01, 0.01, -0.01, 0.01, -0.01],
        "normal_velocity": [3.0, -2.5, -0.01, 0.01, -0.01, 0.01, -0.01, 0.01],
        "blocked": [0] * 8,
        "resolved_shear": [0.0] * 8,
        "free_volume_deficit": [0.0] * 8,
    })
    values.to_csv(tmp_path / "boundary_tracks.csv", index=False)
    metrics = _boundary_metrics(tmp_path)
    assert np.isclose(metrics["raw_reverse_motion_fraction"], 7.0 / 8.0)
    assert np.isclose(metrics["reverse_motion_fraction"], 0.5)


def test_spatial_motion_correlation_and_burst_ccdf_are_reported():
    boundaries = pd.DataFrame({
        "time": [0.0] * 4,
        "grain_i": [0, 0, 2, 3], "grain_j": [1, 2, 3, 4],
        "normal_velocity": [4.0, 4.0, 1.0, 1.0],
    })
    assert _spatial_motion_correlation(boundaries) > 0
    tracks = pd.DataFrame({
        "grain_id": [0, 0, 0, 1, 1, 1],
        "time": [0.0, 1.0, 2.0] * 2,
        "area": [1.0, 2.0, 5.0, 4.0, 4.0, 6.0],
    })
    ccdf = _burst_size_ccdf([tracks])
    assert ccdf["samples"] == 3
    assert ccdf["probability"][0] == 1.0
    assert np.all(np.diff(ccdf["probability"]) <= 0)
    distributions = _trajectory_distributions([tracks])
    assert distributions["absolute_area_rate"]["samples"] == 4
    assert distributions["burst_area_increment"]["samples"] > 0


def test_event_rate_observation_retains_zero_event_exposure(tmp_path):
    pd.DataFrame({
        "time": [0.0, 0.0, 1.0, 1.0, 1.0],
    }).to_csv(tmp_path / "boundary_tracks.csv", index=False)
    pd.DataFrame({"time": []}).to_csv(tmp_path / "events.csv", index=False)
    count, exposure = _event_rate_observation(tmp_path)
    assert count == 0
    assert np.isclose(exposure, 2.5)

    pd.DataFrame({"time": [0.2, 0.8]}).to_csv(tmp_path / "events.csv", index=False)
    count, exposure = _event_rate_observation(tmp_path)
    assert count == 2
    assert np.isclose(exposure, 2.5)

    (tmp_path / "events.csv").unlink()
    pd.DataFrame({
        "time": [0.2, 0.8],
        "event_type": ["activation_hit", "compatibility_release"],
    }).to_csv(tmp_path / "events.csv.gz", index=False, compression="gzip")
    count, exposure = _event_rate_observation(tmp_path)
    assert count == 1
    assert np.isclose(exposure, 2.5)
    (tmp_path / "events.csv.gz").unlink()

    pd.DataFrame({
        "time": [0.2, 0.2, 0.8, 0.8],
        "event_type": ["activation_hit", "compatibility_release",
                       "climb_exchange", "climb_quota_completion"],
    }).to_csv(tmp_path / "events.csv", index=False)
    count, exposure = _event_rate_observation(tmp_path)
    assert count == 2
    assert np.isclose(exposure, 2.5)


def test_event_diagnostics_separates_primitive_rows_and_climb_resistance(tmp_path):
    pd.DataFrame({
        "time": [1.0, 1.0, 1.0, 2.0, 3.0, 3.0],
        "entity_id": ["gb", "gb", "tj:1-2-3", "gb", "gb", "gb"],
        "event_type": ["activation_hit", "compatibility_release",
                       "tj_compatibility_failure", "climb_nucleation",
                       "climb_exchange", "climb_quota_completion"],
        "instantaneous_rate": [4.0, 4.0, 4.0, 2.0, 1.0, 1.0],
        "shear_strain_increment": [0.0, 0.2, 0.0, 0.0, 0.0, 0.1],
        "volumetric_strain_increment": [0.0, 0.0, 0.0, 0.0, 0.0, 0.3],
        "barrier_type": ["easy", "easy", "easy", "climb", "climb", "climb"],
        "DeltaG0": [0.2, 0.2, 0.2, 0.6, 0.6, 0.6],
        "effective_DeltaG": [0.18, 0.18, 0.27, 0.6, 0.6, 0.6],
        "burgers_vector_b": ["[0.1, 0.0]"] * 6,
    }).to_csv(tmp_path / "events.csv", index=False)
    detail = _event_diagnostics([tmp_path])
    assert detail["primitive_event_counts"] == {
        "activation_hit": 1, "climb_exchange": 1, "climb_nucleation": 1,
    }
    assert detail["release_summary_counts"] == {
        "climb_quota_completion": 1, "compatibility_release": 1,
    }
    assert np.isclose(
        detail["climb_expected_resistance_fraction"]["climb_exchange"], 2 / 3
    )
    residence = detail["event_conditioned_expected_residence"]
    assert np.isclose(residence["activation_hit"]["quantiles"]["q50"], 0.25)
    assert np.isclose(residence["climb_nucleation"]["quantiles"]["q50"], 0.5)
    assert np.isclose(residence["climb_exchange"]["quantiles"]["q50"], 1.0)
    fractions = detail["event_conditioned_resistance_fraction"]
    assert np.isclose(sum(fractions.values()), 1.0)
    assert fractions["climb_exchange"] > fractions["climb_nucleation"]
    assert np.isclose(detail["accumulated_event_strain"]["signed_shear"], 0.3)
    tj = detail["tj_compatibility_failures"]
    assert tj["endpoint_failure_rows"] == 1
    assert tj["completed_gb_mode_events"] == 1
    assert tj["low_barrier_failure_rows"] == 1
    assert tj["unique_tj_entities"] == 1
    assert np.isclose(tj["endpoint_failure_incidence_per_mode_event"], 1.0)
    assert np.isclose(
        tj["residual_energy_barrier_shift_ev"]["quantiles"]["q50"], 0.07
    )
    assert np.isclose(tj["packet_burgers_magnitude"]["quantiles"]["q50"], 0.1)

def test_analytical_limits():
    t = np.arange(4.0)
    assert np.allclose(intrinsic_radius(t, 2, 0.5) ** 2, 4 + t)
    assert np.isclose(poisson_activity(1, 2), 1 - np.exp(-2))
    assert np.isclose(series_activity(0.5, 0.25), 1 / 6)
    assert asymptotic_exponent(1, 3) == 5


def test_additive_mechanistic_growth_fits_recover_class_b_and_c():
    time = np.linspace(0.0, 100.0, 101)
    class_b_radius = crossover_radius_prediction(time, 3.0, 0.8, 0.04, 2.0)
    class_b = fit_crossover_growth(time, class_b_radius)
    assert np.isclose(class_b.intrinsic_constant, 0.8, rtol=2e-3)
    assert np.isclose(class_b.crossover_strength, 0.04, rtol=2e-3)
    assert np.isclose(class_b.size_exponent, 2.0, rtol=2e-3)
    assert class_b.r_squared > 0.999999

    class_c_radius = crossover_radius_prediction(time, 3.0, 0.6, 0.1, 1.0)
    class_c = fit_crossover_growth(time, class_c_radius, size_exponent=1.0)
    assert np.isclose(class_c.intrinsic_constant, 0.6, rtol=2e-3)
    assert np.isclose(class_c.crossover_strength, 0.1, rtol=2e-3)
    assert class_c.r_squared > 0.999999


def test_event_ledger_schema(tmp_path):
    target = tmp_path / "events.csv"
    with EventLedger(target) as ledger:
        ledger.write({"run_id": "x", "event_id": 1, "normal_step_h": 0.2})
    with target.open() as handle:
        row = next(csv.DictReader(handle))
    assert tuple(row) == EVENT_FIELDS
    assert row["normal_step_h"] == "0.2"


def test_gzip_event_ledger_truncates_to_checkpoint_member(tmp_path):
    target = tmp_path / "events.csv.gz"
    ledger = EventLedger(target)
    ledger.write({"run_id": "kept", "event_id": 1})
    offset = ledger.checkpoint()
    ledger.write({"run_id": "discarded", "event_id": 2})
    ledger.close()

    resumed = EventLedger(target)
    resumed.truncate(offset)
    resumed.write({"run_id": "resumed", "event_id": 3})
    resumed.close()
    rows = pd.read_csv(target)
    assert rows["run_id"].tolist() == ["kept", "resumed"]


def test_parquet_event_ledger_truncates_to_checkpoint_part(tmp_path):
    target = tmp_path / "events.parquet"
    ledger = EventLedger(target)
    ledger.write({"run_id": "kept", "event_id": "gb:1-2:0:1", "step": 1})
    offset = ledger.checkpoint()
    ledger.write({"run_id": "discarded", "event_id": "gb:1-2:0:2", "step": 2})
    ledger.checkpoint()
    ledger.close()

    resumed = EventLedger(target)
    resumed.truncate(offset)
    resumed.write({
        "run_id": "resumed", "event_id": "gb:1-2:0:3", "step": 2,
        "position": np.asarray([1.5, 2.5]),
    })
    resumed.checkpoint()
    resumed.close()

    assert event_ledger_path(tmp_path) == target
    assert event_ledger_has_rows(target)
    rows = read_event_ledger(target)
    assert rows["run_id"].tolist() == ["kept", "resumed"]
    assert rows["step"].tolist() == [1, 2]
    assert rows["position"].iloc[1] == "[1.5, 2.5]"
    projected = read_event_ledger(target, columns=["event_type", "time", "step"])
    assert list(projected.columns) == ["event_type", "time", "step"]
    assert projected["step"].tolist() == [1, 2]


def test_manifest_can_pin_launch_revision(tmp_path):
    target = tmp_path / "manifest.json"
    write_manifest(target, {"regime": "test"}, "running", code_sha="launch-sha")
    write_manifest(target, {"regime": "test"}, "completed", code_sha="launch-sha")
    manifest = json.loads(target.read_text())
    assert manifest["git_sha"] == "launch-sha"
    assert manifest["status"] == "completed"
