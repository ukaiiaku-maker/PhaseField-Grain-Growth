from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from grain_growth_pf.analysis.activation_energy import fit_activation_energy
from grain_growth_pf.analysis.analytical_models import fit_crossover_growth
from grain_growth_pf.analysis.grain_tracks import ensemble_radius, load_tracks
from grain_growth_pf.analysis.growth_law import (
    fit_common_exponent,
    fit_growth_law,
    fit_growth_law_fixed_exponent,
    scan_growth_exponent,
)
from grain_growth_pf.analysis.jerkiness import jerkiness_metrics
from grain_growth_pf.io.provenance import git_sha


SUMMARY_COLUMNS = [
    "regime", "temperature", "n", "n_ci_low", "n_ci_high", "K", "K_ci",
    "Q_app", "Q_app_ci", "jerkiness_CV", "Fano", "burstiness",
    "reverse_motion_fraction", "velocity_curvature_R2", "pinned_fraction",
    "number_of_events", "number_of_realizations", "Git_SHA",
]

RADIUS_MEASURES = ("R_A", "R_mean", "R_median", "R_rms", "R_perimeter")

CLASS_B_REGIMES = {
    "G1", "G2", "G3", "T1", "T2", "T3", "C1", "C2",
    "P1", "P2", "P3", "P4", "E0", "E1", "E2", "J1", "J2", "J3",
}
CLASS_C_REGIMES = {"C3", "C4", "C5"}
CLASS_D_REGIMES = {"S2", "S3", "SC1", "SC2", "SC3", "SC4", "P5"}


def _growth_window_arrays(run_dirs: list[Path], measure: str = "R_A") -> tuple[np.ndarray, np.ndarray, dict]:
    """Return aligned per-realization radii in the topology-selected window."""
    if measure not in RADIUS_MEASURES:
        raise ValueError(f"unknown radius measure {measure}")
    aligned = None
    for index, path in enumerate(run_dirs):
        _, _, radius = _run_observables(path)
        renamed = radius[["time", measure, "grain_count"]].rename(columns={
            measure: f"R_{index}", "grain_count": f"N_{index}",
        })
        aligned = renamed if aligned is None else aligned.merge(renamed, on="time", how="inner")
    assert aligned is not None
    aligned = aligned.sort_values("time").reset_index(drop=True)
    count_columns = [f"N_{index}" for index in range(len(run_dirs))]
    mean_count = aligned[count_columns].to_numpy(float).mean(axis=1)
    start, end, reason = _fit_window(mean_count)
    time = aligned["time"].to_numpy(float)[start:end]
    radii = aligned[[f"R_{index}" for index in range(len(run_dirs))]].to_numpy(float).T[:, start:end]
    metadata = {
        "selection": reason,
        "initial_mean_grain_count": float(mean_count[0]),
        "end_mean_grain_count": float(mean_count[end - 1]),
        "samples": len(time),
    }
    return time, radii, metadata


def _run_observables(run_dir: Path) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    manifest = json.loads((run_dir / "manifest.json").read_text())
    per_grain = load_tracks(run_dir / "grain_tracks.csv")
    radius = ensemble_radius(per_grain)[
        ["time", "step", *RADIUS_MEASURES, "grain_count"]
    ]
    return manifest, per_grain, radius


def _fit_window(mean_count: np.ndarray) -> tuple[int, int, str]:
    """Choose a recorded topology-based post-transient/pre-finite-size window."""
    initial = float(mean_count[0])
    start_hits = np.flatnonzero(mean_count <= 0.95 * initial)
    start = int(start_hits[0]) if len(start_hits) else max(1, int(0.1 * len(mean_count)))
    end_hits = np.flatnonzero(mean_count < max(20.0, 0.30 * initial))
    end = int(end_hits[0]) if len(end_hits) else len(mean_count)
    if end - start < 8:
        start, end = max(1, int(0.1 * len(mean_count))), len(mean_count)
        reason = "fallback_10pct_to_available_end"
    elif not len(end_hits):
        reason = "five_pct_loss_to_available_end"
    else:
        reason = "five_pct_loss_to_thirty_pct_population"
    return start, end, reason


def _trajectory_metrics(per_grain: pd.DataFrame) -> dict[str, float]:
    metrics = []
    for _, grain in per_grain.groupby("grain_id"):
        grain = grain.sort_values("time")
        if len(grain) < 3:
            continue
        time = grain["time"].to_numpy(float)
        area = grain["area"].to_numpy(float)
        metrics.append(jerkiness_metrics(time, area))
    if not metrics:
        return {}
    keys = set().union(*(metric.keys() for metric in metrics))
    return {
        key: float(np.mean([metric[key] for metric in metrics if key in metric]))
        for key in keys
    }


def _boundary_metrics(run_dir: Path) -> dict[str, float]:
    path = run_dir / "boundary_tracks.csv"
    if not path.exists() or path.stat().st_size == 0:
        return {}
    boundaries = pd.read_csv(path)
    if boundaries.empty:
        return {}
    curvature = boundaries["curvature"].to_numpy(float)
    velocity = boundaries["normal_velocity"].to_numpy(float)
    valid = np.isfinite(curvature) & np.isfinite(velocity) & (np.abs(curvature) > 1e-12) & (np.abs(velocity) > 1e-12)
    reverse = float(np.mean(curvature[valid] * velocity[valid] < 0)) if np.any(valid) else np.nan
    if np.count_nonzero(valid) > 2:
        correlation = np.corrcoef(curvature[valid], velocity[valid])[0, 1]
        curvature_r2 = float(correlation**2) if np.isfinite(correlation) else np.nan
    else:
        curvature_r2 = np.nan
    pinned = float(np.mean(boundaries["blocked"].to_numpy(float)))
    result = {
        "reverse_motion_fraction": reverse,
        "velocity_curvature_R2": curvature_r2,
        "pinned_fraction": pinned,
    }
    for column, name in (
        ("resolved_shear", "velocity_internal_shear_correlation"),
        ("free_volume_deficit", "velocity_free_volume_correlation"),
    ):
        state = boundaries[column].to_numpy(float)
        state_valid = np.isfinite(state) & np.isfinite(velocity)
        if np.count_nonzero(state_valid) > 2 and np.std(state[state_valid]) and np.std(velocity[state_valid]):
            result[name] = float(np.corrcoef(velocity[state_valid], state[state_valid])[0, 1])
        else:
            result[name] = np.nan
    return result


def _neighbor_growth_correlation(per_grain: pd.DataFrame) -> float:
    rates, neighbors = [], []
    for _, grain in per_grain.groupby("grain_id"):
        grain = grain.sort_values("time")
        if len(grain) < 2:
            continue
        elapsed = np.diff(grain["time"].to_numpy(float))
        valid = elapsed > 0
        rates.extend((np.diff(grain["area"].to_numpy(float))[valid] / elapsed[valid]).tolist())
        neighbors.extend(grain["neighbors"].to_numpy(float)[:-1][valid].tolist())
    if len(rates) < 3 or not np.std(rates) or not np.std(neighbors):
        return np.nan
    return float(np.corrcoef(rates, neighbors)[0, 1])


def _nanmean_metric(metrics: list[dict[str, float]], key: str) -> float:
    values = [metric.get(key, np.nan) for metric in metrics]
    return float(np.nanmean(values)) if any(np.isfinite(values)) else np.nan


def _nanmean_values(values: list[float]) -> float:
    return float(np.nanmean(values)) if any(np.isfinite(values)) else np.nan


def _event_statistics(run_dir: Path) -> tuple[int, float]:
    path = run_dir / "events.csv"
    if not path.exists() or path.stat().st_size == 0:
        return 0, np.nan
    events = pd.read_csv(path)
    if events.empty:
        return 0, np.nan
    if "time" not in events:
        return len(events), np.nan
    tracks = pd.read_csv(run_dir / "grain_tracks.csv", usecols=["time"])
    duration = int(np.ceil(tracks["time"].max())) + 1
    counts = np.bincount(np.floor(events["time"].to_numpy(float)).astype(int), minlength=duration)
    fano = float(np.var(counts) / np.mean(counts)) if np.mean(counts) else np.nan
    return len(events), fano


def analyze_group(run_dirs: list[Path], bootstrap_samples: int = 500) -> tuple[dict[str, object], dict]:
    loaded = [_run_observables(path) for path in run_dirs]
    manifests = [item[0] for item in loaded]
    config = manifests[0]["config"]
    aligned = None
    for index, (_, _, radius) in enumerate(loaded):
        renamed = radius.rename(columns={
            **{measure: f"{measure}_{index}" for measure in RADIUS_MEASURES},
            "grain_count": f"N_{index}",
        })
        renamed = renamed.drop(columns="step")
        aligned = renamed if aligned is None else aligned.merge(renamed, on="time", how="inner")
    assert aligned is not None
    aligned = aligned.sort_values("time").reset_index(drop=True)
    radius_columns = [f"R_A_{i}" for i in range(len(loaded))]
    count_columns = [f"N_{i}" for i in range(len(loaded))]
    radii = aligned[radius_columns].to_numpy(float).T
    mean_count = aligned[count_columns].to_numpy(float).mean(axis=1)
    start, end, window_reason = _fit_window(mean_count)
    time = aligned["time"].to_numpy(float)[start:end]
    fit_radii = radii[:, start:end]
    fit = fit_growth_law(time, fit_radii.mean(axis=0), transient_fraction=0.0)
    parabolic_fit = fit_growth_law_fixed_exponent(
        time, fit_radii.mean(axis=0), 2.0, transient_fraction=0.0
    )
    profile = scan_growth_exponent(time, fit_radii.mean(axis=0))

    digest = hashlib.sha256(f"{config['regime']}:{config['pf']['temperature']}".encode()).digest()
    rng = np.random.default_rng(int.from_bytes(digest[:8], "little"))
    bootstrap_n, bootstrap_k = [], []
    for _ in range(bootstrap_samples):
        selection = rng.integers(0, len(fit_radii), len(fit_radii))
        sample_fit = fit_growth_law(time, fit_radii[selection].mean(axis=0), transient_fraction=0.0)
        bootstrap_n.append(sample_fit.exponent)
        # K has exponent-dependent units. Its interval is therefore evaluated
        # at the ensemble best-fit n while n itself is bootstrapped freely.
        sample_radius = fit_radii[selection].mean(axis=0)
        bootstrap_k.append(fit_growth_law_fixed_exponent(
            time, sample_radius, fit.exponent, transient_fraction=0.0
        ).coefficient)
    n_low, n_high = np.quantile(bootstrap_n, [0.025, 0.975])
    k_low, k_high = np.quantile(bootstrap_k, [0.025, 0.975])

    trajectory_metrics = [_trajectory_metrics(item[1]) for item in loaded]
    boundary_metrics = [_boundary_metrics(path) for path in run_dirs]
    event_metrics = [_event_statistics(path) for path in run_dirs]
    profile_best = int(np.argmin(profile.normalized_rmse))
    at_bound = bool(fit.exponent <= 1.01 or fit.exponent >= 5.99)
    row = {
        "regime": config["regime"], "temperature": config["pf"]["temperature"],
        "n": fit.exponent, "n_ci_low": float(n_low), "n_ci_high": float(n_high),
        "K": fit.coefficient, "K_ci": float((k_high - k_low) / 2),
        "Q_app": np.nan, "Q_app_ci": np.nan,
        "jerkiness_CV": _nanmean_metric(trajectory_metrics, "jerkiness_CV"),
        "Fano": float(np.nanmean([m[1] for m in event_metrics])) if any(np.isfinite(m[1]) for m in event_metrics) else np.nan,
        "burstiness": _nanmean_metric(trajectory_metrics, "burstiness"),
        "reverse_motion_fraction": _nanmean_metric(boundary_metrics, "reverse_motion_fraction"),
        "velocity_curvature_R2": _nanmean_metric(boundary_metrics, "velocity_curvature_R2"),
        "pinned_fraction": _nanmean_metric(boundary_metrics, "pinned_fraction"),
        "number_of_events": int(sum(m[0] for m in event_metrics)),
        "number_of_realizations": len(run_dirs), "Git_SHA": manifests[0]["git_sha"],
    }
    diagnostics = {
        "regime": config["regime"], "temperature": config["pf"]["temperature"],
        "fit_window": {"start": fit.fit_start, "end": fit.fit_end, "selection": window_reason,
                       "samples": len(time), "initial_mean_grain_count": float(mean_count[0]),
                       "end_mean_grain_count": float(mean_count[end - 1])},
        "ensemble_fit": {"n": fit.exponent, "K": fit.coefficient, "intercept": fit.intercept,
                         "r_squared": fit.r_squared,
                         "residual_autocorrelation": fit.residual_autocorrelation,
                         "at_search_bound": at_bound},
        "bootstrap": {"samples": bootstrap_samples, "n_95pct": [float(n_low), float(n_high)],
                      "K_95pct": [float(k_low), float(k_high)]},
        "profile_scan": {"best_n": float(profile.exponents[profile_best]),
                         "best_normalized_rmse": float(profile.normalized_rmse[profile_best]),
                         "n_grid": profile.exponents.tolist(),
                         "normalized_rmse": profile.normalized_rmse.tolist(),
                         "residual_autocorrelation": profile.residual_autocorrelation.tolist()},
        "mechanistic_comparators": {
            "intrinsic_parabolic_n2": {
                "K": parabolic_fit.coefficient,
                "intercept": parabolic_fit.intercept,
                "r_squared": parabolic_fit.r_squared,
                "residual_autocorrelation": parabolic_fit.residual_autocorrelation,
                "delta_r_squared_vs_free_n": parabolic_fit.r_squared - fit.r_squared,
            }
        },
        "radius_measure_fits": {},
        "population_band_sensitivity": [],
        "intermittency": {
            key: _nanmean_metric(trajectory_metrics, key)
            for key in (
                "jerkiness_CV", "stationary_fraction", "motion_top_1pct",
                "motion_top_5pct", "motion_top_10pct", "burstiness",
            )
        },
        "correlations": {
            "velocity_internal_shear": _nanmean_metric(
                boundary_metrics, "velocity_internal_shear_correlation"
            ),
            "velocity_free_volume_deficit": _nanmean_metric(
                boundary_metrics, "velocity_free_volume_correlation"
            ),
            "neighbor_number_growth_rate": _nanmean_values([
                _neighbor_growth_correlation(item[1]) for item in loaded
            ]),
        },
    }
    regime = str(config["regime"])
    theory_class = (
        "class_b" if regime in CLASS_B_REGIMES else
        "class_c" if regime in CLASS_C_REGIMES else
        "class_d" if regime in CLASS_D_REGIMES else
        "intrinsic_or_unclassified"
    )
    diagnostics["mechanistic_comparators"]["source_theory_class"] = theory_class
    if theory_class in {"class_b", "class_d"}:
        class_b = fit_crossover_growth(time, fit_radii.mean(axis=0))
        diagnostics["mechanistic_comparators"]["class_b_additive"] = {
            "intrinsic_K": class_b.intrinsic_constant,
            "crossover_strength": class_b.crossover_strength,
            "hazard_size_exponent": class_b.size_exponent,
            "asymptotic_growth_exponent": class_b.size_exponent + 2.0,
            "r_squared": class_b.r_squared,
            "normalized_rmse": class_b.normalized_rmse,
            "parameter_at_bound": class_b.parameter_at_bound,
        }
    if theory_class in {"class_c", "class_d"}:
        class_c = fit_crossover_growth(time, fit_radii.mean(axis=0), size_exponent=1.0)
        diagnostics["mechanistic_comparators"]["class_c_exchange"] = {
            "intrinsic_K": class_c.intrinsic_constant,
            "crossover_strength": class_c.crossover_strength,
            "crossover_radius": 1.0 / class_c.crossover_strength,
            "r_squared": class_c.r_squared,
            "normalized_rmse": class_c.normalized_rmse,
            "parameter_at_bound": class_c.parameter_at_bound,
        }
    for measure in RADIUS_MEASURES:
        columns = [f"{measure}_{index}" for index in range(len(loaded))]
        measure_radius = aligned[columns].to_numpy(float).T[:, start:end].mean(axis=0)
        measure_fit = fit_growth_law(time, measure_radius, transient_fraction=0.0)
        diagnostics["radius_measure_fits"][measure] = {
            "n": measure_fit.exponent,
            "K": measure_fit.coefficient,
            "r_squared": measure_fit.r_squared,
            "residual_autocorrelation": measure_fit.residual_autocorrelation,
        }
    full_time = aligned["time"].to_numpy(float)
    full_radius = radii.mean(axis=0)
    initial_count = float(mean_count[0])
    for upper, lower in ((0.95, 0.75), (0.75, 0.60), (0.60, 0.50), (0.50, 0.40), (0.40, 0.30)):
        selection = (mean_count <= upper * initial_count) & (mean_count >= lower * initial_count)
        if np.count_nonzero(selection) < 8:
            continue
        band_fit = fit_growth_law(
            full_time[selection], full_radius[selection], transient_fraction=0.0
        )
        diagnostics["population_band_sensitivity"].append({
            "upper_fraction": upper,
            "lower_fraction": lower,
            "samples": int(np.count_nonzero(selection)),
            "n": band_fit.exponent,
            "K": band_fit.coefficient,
            "r_squared": band_fit.r_squared,
            "residual_autocorrelation": band_fit.residual_autocorrelation,
        })
    return row, diagnostics


def analyze_run(run_dir: str | Path) -> dict[str, object]:
    """Analyze one run using the same estimator as a one-member ensemble."""
    return analyze_group([Path(run_dir)], bootstrap_samples=1)[0]


def analyze_campaign(campaign_dir: str | Path, output: str | Path | None = None,
                     bootstrap_samples: int = 500) -> pd.DataFrame:
    campaign_dir = Path(campaign_dir)
    manifest = json.loads((campaign_dir / "campaign_manifest.json").read_text())
    grouped: dict[tuple[str, float], list[Path]] = {}
    for raw_path in manifest["runs"]:
        path = Path(raw_path)
        run_manifest = json.loads((path / "manifest.json").read_text())
        config = run_manifest["config"]
        key = (config["regime"], float(config["pf"]["temperature"]))
        grouped.setdefault(key, []).append(path)

    rows, diagnostics = [], []
    for paths in grouped.values():
        row, detail = analyze_group(paths, bootstrap_samples=bootstrap_samples)
        run_manifests = [json.loads((path / "manifest.json").read_text()) for path in paths]
        detail["provenance"] = {
            "campaign": str(campaign_dir),
            "analysis_git_sha": git_sha(),
            "simulation_git_shas": sorted({item["git_sha"] for item in run_manifests}),
            "config_sha256s": sorted({
                item["config_sha256"] for item in run_manifests if item.get("config_sha256")
            }),
        }
        rows.append(row)
        diagnostics.append(detail)
    summary = pd.DataFrame(rows, columns=SUMMARY_COLUMNS)
    detail_by_key = {
        (detail["regime"], float(detail["temperature"])): detail for detail in diagnostics
    }
    for regime, indices in summary.groupby("regime").groups.items():
        subset = summary.loc[indices]
        if len(subset) < 4:
            continue
        temperatures = np.sort(subset["temperature"].to_numpy(float))
        series_times, series_radii, window_metadata = [], [], []
        for temperature in temperatures:
            time, radii, metadata = _growth_window_arrays(grouped[(regime, temperature)])
            series_times.append(time)
            series_radii.append(radii)
            window_metadata.append(metadata)
        common = fit_common_exponent(
            series_times, [radii.mean(axis=0) for radii in series_radii]
        )
        digest = hashlib.sha256(f"temperature-series:{regime}".encode()).digest()
        rng = np.random.default_rng(int.from_bytes(digest[:8], "little"))
        bootstrap_n, bootstrap_k, bootstrap_q = [], [], []
        for _ in range(bootstrap_samples):
            sampled_mean_radii = []
            for radii in series_radii:
                selection = rng.integers(0, len(radii), len(radii))
                sampled_mean_radii.append(radii[selection].mean(axis=0))
            sample_fit = fit_common_exponent(series_times, sampled_mean_radii)
            bootstrap_n.append(sample_fit.exponent)
            bootstrap_k.append(sample_fit.coefficients)
            if np.all(sample_fit.coefficients > 0):
                bootstrap_q.append(fit_activation_energy(
                    temperatures, sample_fit.coefficients
                ).activation_energy_ev)
        bootstrap_k_array = np.asarray(bootstrap_k)
        n_low, n_high = np.quantile(bootstrap_n, [0.025, 0.975])
        activation = fit_activation_energy(temperatures, common.coefficients)
        q_low, q_high = np.quantile(bootstrap_q, [0.025, 0.975])
        for position, temperature in enumerate(temperatures):
            row_index = subset.index[np.isclose(subset["temperature"], temperature)][0]
            k_low, k_high = np.quantile(bootstrap_k_array[:, position], [0.025, 0.975])
            summary.loc[row_index, ["n", "n_ci_low", "n_ci_high"]] = [
                common.exponent, n_low, n_high
            ]
            summary.loc[row_index, ["K", "K_ci"]] = [
                common.coefficients[position], (k_high - k_low) / 2.0
            ]
            summary.loc[row_index, ["Q_app", "Q_app_ci"]] = [
                activation.activation_energy_ev, (q_high - q_low) / 2.0
            ]
            detail_by_key[(regime, float(temperature))]["temperature_series_fit"] = {
                "common_n": common.exponent,
                "common_n_95pct": [float(n_low), float(n_high)],
                "normalized_rmse": common.normalized_rmse,
                "temperatures": temperatures.tolist(),
                "coefficients": common.coefficients.tolist(),
                "activation_energy_ev": activation.activation_energy_ev,
                "activation_energy_95pct": [float(q_low), float(q_high)],
                "bootstrap_samples": bootstrap_samples,
                "window_by_temperature": window_metadata,
            }
    target = Path(output) if output else campaign_dir / "mechanism_summary.csv"
    summary.to_csv(target, index=False)
    target.with_name(f"{target.stem}_diagnostics.json").write_text(json.dumps(diagnostics, indent=2))
    return summary
