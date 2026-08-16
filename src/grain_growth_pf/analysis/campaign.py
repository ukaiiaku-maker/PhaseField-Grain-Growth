from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from grain_growth_pf.analysis.activation_energy import (
    fit_activation_energy,
    local_activation_energies,
)
from grain_growth_pf.analysis.analytical_models import fit_crossover_growth
from grain_growth_pf.analysis.grain_tracks import ensemble_radius, load_tracks
from grain_growth_pf.analysis.growth_law import (
    fit_common_exponent,
    fit_growth_law,
    fit_growth_law_fixed_exponent,
    scan_growth_exponent,
)
from grain_growth_pf.analysis.jerkiness import jerkiness_metrics
from grain_growth_pf.io.event_ledger import event_ledger_path
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
ACTIVATION_EVENT_TYPES = {
    "activation_hit", "tj_activation_hit", "climb_nucleation",
    "climb_exchange", "climb_transport",
}


def _activation_rows(events: pd.DataFrame) -> pd.DataFrame:
    """Prefer primitive stochastic passages over duplicate release summaries."""
    if "event_type" not in events:
        return events
    selected = events[events["event_type"].isin(ACTIVATION_EVENT_TYPES)]
    return selected if not selected.empty else events


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
        "relative_radius_change": float(
            (radii[:, -1].mean() - radii[:, 0].mean())
            / max(radii[:, 0].mean(), np.finfo(float).tiny)
        ),
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
    raw_reverse = float(np.mean(curvature[valid] * velocity[valid] < 0)) if np.any(valid) else np.nan
    if np.any(valid):
        curvature_threshold = float(np.quantile(np.abs(curvature[valid]), 0.75))
        velocity_threshold = float(np.quantile(np.abs(velocity[valid]), 0.75))
        active = valid & (np.abs(curvature) >= curvature_threshold)
        active &= np.abs(velocity) >= velocity_threshold
        reverse = float(np.mean(curvature[active] * velocity[active] < 0)) if np.any(active) else np.nan
    else:
        active = np.zeros_like(valid)
        reverse = np.nan
    if np.count_nonzero(valid) > 2:
        correlation = np.corrcoef(curvature[valid], velocity[valid])[0, 1]
        curvature_r2 = float(correlation**2) if np.isfinite(correlation) else np.nan
    else:
        curvature_r2 = np.nan
    pinned = float(np.mean(boundaries["blocked"].to_numpy(float)))
    result = {
        "reverse_motion_fraction": reverse,
        "raw_reverse_motion_fraction": raw_reverse,
        "active_boundary_fraction": float(np.mean(active)),
        "velocity_curvature_R2": curvature_r2,
        "pinned_fraction": pinned,
        "simultaneous_motion_spatial_correlation": _spatial_motion_correlation(boundaries),
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


def _spatial_motion_correlation(boundaries: pd.DataFrame,
                                maximum_time_samples: int = 200) -> float:
    """Topological Moran-like correlation of simultaneous GB speed.

    Boundary domains are spatial neighbors when they share a grain. The metric
    compares centered speed products for neighboring domains with the global
    speed variance at the same recorded time.
    """
    required = {"time", "grain_i", "grain_j", "normal_velocity"}
    if not required.issubset(boundaries.columns) or boundaries.empty:
        return np.nan
    grouped = boundaries.groupby("time", sort=False)
    times = np.asarray(list(grouped.indices), dtype=float)
    if len(times) > maximum_time_samples:
        times = times[np.unique(np.linspace(
            0, len(times) - 1, maximum_time_samples, dtype=int
        ))]
    correlations = []
    for time_value in times:
        frame = grouped.get_group(time_value)
        velocity = np.abs(frame["normal_velocity"].to_numpy(float))
        finite = np.isfinite(velocity)
        if np.count_nonzero(finite) < 3 or np.var(velocity[finite]) <= 0:
            continue
        frame = frame.iloc[np.flatnonzero(finite)]
        centered = velocity[finite] - np.mean(velocity[finite])
        grain_ids = np.concatenate((
            frame["grain_i"].to_numpy(int), frame["grain_j"].to_numpy(int)
        ))
        values = np.concatenate((centered, centered))
        order = np.argsort(grain_ids, kind="stable")
        sorted_ids, sorted_values = grain_ids[order], values[order]
        starts = np.r_[0, np.flatnonzero(np.diff(sorted_ids)) + 1]
        counts = np.diff(np.r_[starts, len(sorted_ids)])
        sums = np.add.reduceat(sorted_values, starts)
        sums_squared = np.add.reduceat(sorted_values**2, starts)
        pair_counts = counts * (counts - 1) / 2.0
        valid_groups = pair_counts > 0
        total_pairs = float(pair_counts[valid_groups].sum())
        if total_pairs:
            pair_products = 0.5 * (sums**2 - sums_squared)
            correlations.append(float(
                pair_products[valid_groups].sum()
                / (total_pairs * np.var(centered))
            ))
    return float(np.mean(correlations)) if correlations else np.nan


def _burst_size_ccdf(per_grain_frames: list[pd.DataFrame]) -> dict[str, object]:
    sizes = []
    for frame in per_grain_frames:
        for _, grain in frame.groupby("grain_id"):
            increments = np.abs(np.diff(grain.sort_values("time")["area"].to_numpy(float)))
            sizes.extend(increments[increments > 1e-12].tolist())
    if not sizes:
        return {"samples": 0, "size": [], "probability": []}
    values = np.asarray(sizes, dtype=float)
    thresholds = np.unique(np.quantile(values, np.linspace(0.0, 1.0, 65)))
    return {
        "samples": int(len(values)),
        "size": thresholds.tolist(),
        "probability": [float(np.mean(values >= threshold)) for threshold in thresholds],
    }


def _trajectory_distributions(per_grain_frames: list[pd.DataFrame]) -> dict[str, object]:
    """Compact grain-scale rate, waiting, burst-duration, and burst-size data."""
    rates, waits, durations, burst_sizes = [], [], [], []
    for frame in per_grain_frames:
        for _, grain in frame.groupby("grain_id"):
            grain = grain.sort_values("time")
            time = grain["time"].to_numpy(float)
            area = grain["area"].to_numpy(float)
            if len(time) < 3:
                continue
            dt = np.diff(time)
            valid = dt > 0
            if np.count_nonzero(valid) < 2:
                continue
            magnitude = np.abs(np.diff(area)[valid] / dt[valid])
            interval_time = time[1:][valid]
            interval_dt = dt[valid]
            increments = np.abs(np.diff(area)[valid])
            rates.extend(magnitude.tolist())
            positive = magnitude[magnitude > 1e-12]
            if not len(positive):
                continue
            threshold = float(np.quantile(positive, 0.9))
            active = magnitude >= threshold
            starts = np.flatnonzero(active & ~np.r_[False, active[:-1]])
            ends = np.flatnonzero(active & ~np.r_[active[1:], False])
            if len(starts) > 1:
                waits.extend((interval_time[starts[1:]] - interval_time[ends[:-1]]).tolist())
            for start, end in zip(starts, ends):
                durations.append(float(interval_dt[start:end + 1].sum()))
                burst_sizes.append(float(increments[start:end + 1].sum()))
    return {
        "absolute_area_rate": _quantiles(rates),
        "interburst_waiting_time": _quantiles(waits),
        "burst_duration": _quantiles(durations),
        "burst_area_increment": _quantiles(burst_sizes),
        "burst_definition": "contiguous intervals at or above each grain's 90th percentile positive absolute area rate",
    }


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
    path = event_ledger_path(run_dir)
    if not path.exists() or path.stat().st_size == 0:
        return 0, np.nan
    events = pd.read_csv(path)
    if events.empty:
        return 0, np.nan
    events = _activation_rows(events)
    if "time" not in events:
        return len(events), np.nan
    tracks = pd.read_csv(run_dir / "grain_tracks.csv", usecols=["time"])
    duration = int(np.ceil(tracks["time"].max())) + 1
    counts = np.bincount(np.floor(events["time"].to_numpy(float)).astype(int), minlength=duration)
    fano = float(np.var(counts) / np.mean(counts)) if np.mean(counts) else np.nan
    return len(events), fano


def _event_rate_observation(run_dir: Path) -> tuple[int, float]:
    """Return events and integrated GB-domain exposure for censoring-aware rates."""
    event_path = event_ledger_path(run_dir)
    boundary_path = run_dir / "boundary_tracks.csv"
    if not boundary_path.exists() or boundary_path.stat().st_size == 0:
        return 0, 0.0
    times = pd.read_csv(boundary_path, usecols=["time"])["time"]
    counts = times.value_counts(sort=False).sort_index()
    if len(counts) < 2:
        exposure = 0.0
    else:
        exposure = float(np.trapezoid(
            counts.to_numpy(float), counts.index.to_numpy(float)
        ))
    if not event_path.exists() or event_path.stat().st_size == 0:
        return 0, exposure
    events = _activation_rows(pd.read_csv(event_path))
    return len(events), exposure


def _quantiles(values: list[float] | np.ndarray) -> dict[str, object]:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if not len(array):
        return {"samples": 0, "quantiles": {}}
    levels = np.asarray([0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0])
    return {
        "samples": int(len(array)),
        "quantiles": {
            f"q{int(level * 100):02d}": float(value)
            for level, value in zip(levels, np.quantile(array, levels))
        },
    }


def _event_diagnostics(run_dirs: list[Path]) -> dict[str, object]:
    """Summarize primitive first passages without mixing independent clocks."""
    frames = []
    for run_index, run_dir in enumerate(run_dirs):
        path = event_ledger_path(run_dir)
        if not path.exists() or path.stat().st_size == 0:
            continue
        frame = pd.read_csv(path)
        if frame.empty:
            continue
        frame["realization"] = run_index
        frames.append(frame)
    if not frames:
        return {"primitive_event_counts": {}, "waiting_times": _quantiles([])}
    all_events = pd.concat(frames, ignore_index=True)
    primitive = _activation_rows(all_events)
    counts = primitive["event_type"].value_counts().sort_index().to_dict()
    waiting_times: list[float] = []
    for _, entity in primitive.groupby(["realization", "entity_id"], dropna=False):
        differences = np.diff(np.sort(entity["time"].to_numpy(float)))
        waiting_times.extend(differences[differences > 0].tolist())
    stage_residence = {}
    stage_rows = primitive[primitive["event_type"].isin({
        "climb_nucleation", "climb_exchange", "climb_transport",
    })]
    resistance = {}
    for event_type, rows in stage_rows.groupby("event_type"):
        rates = pd.to_numeric(rows["instantaneous_rate"], errors="coerce").to_numpy(float)
        rates = rates[np.isfinite(rates) & (rates > 0)]
        stage_residence[event_type] = _quantiles(1.0 / rates)
        resistance[event_type] = float(np.mean(1.0 / rates)) if len(rates) else np.nan
    finite_resistance = {key: value for key, value in resistance.items() if np.isfinite(value)}
    total_resistance = sum(finite_resistance.values())
    resistance_fraction = {
        key: value / total_resistance for key, value in finite_resistance.items()
    } if total_resistance else {}
    shear_increment = (
        pd.to_numeric(all_events["shear_strain_increment"], errors="coerce")
        if "shear_strain_increment" in all_events else pd.Series(dtype=float)
    )
    volumetric_increment = (
        pd.to_numeric(all_events["volumetric_strain_increment"], errors="coerce")
        if "volumetric_strain_increment" in all_events else pd.Series(dtype=float)
    )
    return {
        "primitive_event_counts": {str(key): int(value) for key, value in counts.items()},
        "waiting_times_by_entity": _quantiles(waiting_times),
        "climb_expected_stage_residence": stage_residence,
        "climb_expected_resistance_fraction": resistance_fraction,
        "release_summary_counts": {
            str(key): int(value) for key, value in
            all_events[~all_events.index.isin(primitive.index)]["event_type"]
            .value_counts().sort_index().items()
        },
        "accumulated_event_strain": {
            "signed_shear": float(shear_increment.fillna(0.0).sum()),
            "signed_volumetric": float(volumetric_increment.fillna(0.0).sum()),
        },
    }


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
    bootstrap_iterations = 1 if len(fit_radii) == 1 else bootstrap_samples
    for _ in range(bootstrap_iterations):
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
        "bootstrap": {"samples": bootstrap_iterations,
                      "requested_samples": bootstrap_samples,
                      "degenerate_single_realization": len(fit_radii) == 1,
                      "n_95pct": [float(n_low), float(n_high)],
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
        "burst_size_ccdf": _burst_size_ccdf([item[1] for item in loaded]),
        "trajectory_distributions": _trajectory_distributions([item[1] for item in loaded]),
        "event_diagnostics": _event_diagnostics(run_dirs),
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
            "raw_reverse_motion_fraction": _nanmean_metric(
                boundary_metrics, "raw_reverse_motion_fraction"
            ),
            "active_boundary_fraction": _nanmean_metric(
                boundary_metrics, "active_boundary_fraction"
            ),
            "simultaneous_boundary_motion_spatial": _nanmean_metric(
                boundary_metrics, "simultaneous_motion_spatial_correlation"
            ),
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
                     bootstrap_samples: int = 500,
                     require_completed: bool = True) -> pd.DataFrame:
    campaign_dir = Path(campaign_dir)
    analysis_sha = git_sha()
    manifest = json.loads((campaign_dir / "campaign_manifest.json").read_text())
    if require_completed and manifest.get("status") != "completed":
        raise ValueError(
            f"campaign is not complete ({manifest.get('status', 'missing status')}): "
            f"{campaign_dir}; pass require_completed=False only for diagnostic analysis"
        )
    grouped: dict[tuple[str, float], list[Path]] = {}
    for raw_path in manifest["runs"]:
        path = Path(raw_path)
        run_manifest_path = path / "manifest.json"
        if not run_manifest_path.exists():
            if require_completed:
                raise ValueError(f"run has not started: {path}")
            continue
        run_manifest = json.loads(run_manifest_path.read_text())
        if require_completed and run_manifest.get("status") != "completed":
            raise ValueError(f"run is not complete: {path}")
        if run_manifest.get("status") != "completed":
            continue
        config = run_manifest["config"]
        key = (config["regime"], float(config["pf"]["temperature"]))
        grouped.setdefault(key, []).append(path)

    rows, diagnostics = [], []
    for paths in grouped.values():
        row, detail = analyze_group(paths, bootstrap_samples=bootstrap_samples)
        run_manifests = [json.loads((path / "manifest.json").read_text()) for path in paths]
        detail["provenance"] = {
            "campaign": str(campaign_dir),
            "analysis_git_sha": analysis_sha,
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
        event_observations = []
        for temperature in temperatures:
            time, radii, metadata = _growth_window_arrays(grouped[(regime, temperature)])
            series_times.append(time)
            series_radii.append(radii)
            window_metadata.append(metadata)
            event_observations.append([
                _event_rate_observation(path) for path in grouped[(regime, temperature)]
            ])
        common = fit_common_exponent(
            series_times, [radii.mean(axis=0) for radii in series_radii]
        )
        minimum_relative_growth = 0.02
        observable = all(
            metadata["relative_radius_change"] >= minimum_relative_growth
            for metadata in window_metadata
        )
        digest = hashlib.sha256(f"temperature-series:{regime}".encode()).digest()
        rng = np.random.default_rng(int.from_bytes(digest[:8], "little"))
        bootstrap_n, bootstrap_k, bootstrap_q, bootstrap_event_q = [], [], [], []
        for _ in range(bootstrap_samples):
            sampled_mean_radii, sampled_event_rates = [], []
            for radii, observations in zip(series_radii, event_observations):
                selection = rng.integers(0, len(radii), len(radii))
                sampled_mean_radii.append(radii[selection].mean(axis=0))
                event_count = sum(observations[index][0] for index in selection)
                exposure = sum(observations[index][1] for index in selection)
                sampled_event_rates.append(event_count / exposure if exposure > 0 else 0.0)
            sample_fit = fit_common_exponent(series_times, sampled_mean_radii)
            bootstrap_n.append(sample_fit.exponent)
            bootstrap_k.append(sample_fit.coefficients)
            if observable and np.all(sample_fit.coefficients > 0):
                bootstrap_q.append(fit_activation_energy(
                    temperatures, sample_fit.coefficients
                ).activation_energy_ev)
            if np.all(np.asarray(sampled_event_rates) > 0):
                bootstrap_event_q.append(fit_activation_energy(
                    temperatures, np.asarray(sampled_event_rates)
                ).activation_energy_ev)
        bootstrap_k_array = np.asarray(bootstrap_k)
        n_low, n_high = np.quantile(bootstrap_n, [0.025, 0.975])
        activation = (
            fit_activation_energy(temperatures, common.coefficients)
            if observable else None
        )
        if activation is not None:
            local_temperature, local_q = local_activation_energies(
                temperatures, common.coefficients
            )
        else:
            local_temperature = local_q = np.asarray([])
        q_interval = (
            np.quantile(bootstrap_q, [0.025, 0.975]).tolist()
            if bootstrap_q else None
        )
        pooled_event_counts = np.asarray([
            sum(observation[0] for observation in observations)
            for observations in event_observations
        ])
        pooled_event_exposure = np.asarray([
            sum(observation[1] for observation in observations)
            for observations in event_observations
        ])
        pooled_event_rates = np.divide(
            pooled_event_counts, pooled_event_exposure,
            out=np.zeros_like(pooled_event_exposure), where=pooled_event_exposure > 0,
        )
        event_activation = (
            fit_activation_energy(temperatures, pooled_event_rates)
            if np.all(pooled_event_rates > 0) else None
        )
        event_q_interval = (
            np.quantile(bootstrap_event_q, [0.025, 0.975]).tolist()
            if bootstrap_event_q else None
        )
        if event_activation is not None:
            event_local_temperature, event_local_q = local_activation_energies(
                temperatures, pooled_event_rates
            )
        else:
            event_local_temperature = event_local_q = np.asarray([])
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
                activation.activation_energy_ev if activation else np.nan,
                (q_interval[1] - q_interval[0]) / 2.0 if q_interval else np.nan,
            ]
            detail_by_key[(regime, float(temperature))]["temperature_series_fit"] = {
                "common_n": common.exponent,
                "common_n_95pct": [float(n_low), float(n_high)],
                "normalized_rmse": common.normalized_rmse,
                "temperatures": temperatures.tolist(),
                "coefficients": common.coefficients.tolist(),
                "activation_energy_ev": (
                    activation.activation_energy_ev if activation else None
                ),
                "activation_energy_95pct": q_interval,
                "arrhenius_r_squared": activation.r_squared if activation else None,
                "arrhenius_standard_error_ev": (
                    activation.standard_error_ev if activation else None
                ),
                "local_activation_midpoint_temperature": local_temperature.tolist(),
                "local_activation_energy_ev": local_q.tolist(),
                "kinetically_observable": observable,
                "minimum_required_relative_radius_change": minimum_relative_growth,
                "bootstrap_samples": bootstrap_samples,
                "window_by_temperature": window_metadata,
                "event_level": {
                    "estimator": "event_count_per_integrated_GB_domain_time",
                    "counts": pooled_event_counts.tolist(),
                    "domain_time_exposure": pooled_event_exposure.tolist(),
                    "rates": pooled_event_rates.tolist(),
                    "activation_energy_ev": (
                        event_activation.activation_energy_ev if event_activation else None
                    ),
                    "activation_energy_95pct": event_q_interval,
                    "arrhenius_r_squared": (
                        event_activation.r_squared if event_activation else None
                    ),
                    "local_activation_midpoint_temperature": (
                        event_local_temperature.tolist()
                    ),
                    "local_activation_energy_ev": event_local_q.tolist(),
                    "bootstrap_samples_with_events_at_all_temperatures": len(
                        bootstrap_event_q
                    ),
                },
            }
    target = Path(output) if output else campaign_dir / "mechanism_summary.csv"
    summary.to_csv(target, index=False)
    target.with_name(f"{target.stem}_diagnostics.json").write_text(json.dumps(diagnostics, indent=2))
    return summary
