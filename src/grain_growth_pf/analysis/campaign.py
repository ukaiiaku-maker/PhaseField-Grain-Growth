from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from grain_growth_pf.analysis.activation_energy import fit_activation_energy
from grain_growth_pf.analysis.grain_tracks import ensemble_radius, load_tracks
from grain_growth_pf.analysis.growth_law import fit_growth_law, scan_growth_exponent
from grain_growth_pf.analysis.jerkiness import jerkiness_metrics


SUMMARY_COLUMNS = [
    "regime", "temperature", "n", "n_ci_low", "n_ci_high", "K", "K_ci",
    "Q_app", "Q_app_ci", "jerkiness_CV", "Fano", "burstiness",
    "reverse_motion_fraction", "velocity_curvature_R2", "pinned_fraction",
    "number_of_events", "number_of_realizations", "Git_SHA",
]


def _run_observables(run_dir: Path) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    manifest = json.loads((run_dir / "manifest.json").read_text())
    per_grain = load_tracks(run_dir / "grain_tracks.csv")
    radius = ensemble_radius(per_grain)[["time", "step", "R_A", "grain_count"]]
    return manifest, per_grain, radius


def _fit_window(mean_count: np.ndarray) -> tuple[int, int, str]:
    """Choose a recorded topology-based post-transient/pre-finite-size window."""
    initial = float(mean_count[0])
    start_hits = np.flatnonzero(mean_count <= 0.95 * initial)
    start = int(start_hits[0]) if len(start_hits) else max(1, int(0.1 * len(mean_count)))
    end_hits = np.flatnonzero(mean_count < max(20.0, 0.60 * initial))
    end = int(end_hits[0]) if len(end_hits) else len(mean_count)
    if end - start < 8:
        start, end = max(1, int(0.1 * len(mean_count))), len(mean_count)
        reason = "fallback_10pct_to_available_end"
    elif not len(end_hits):
        reason = "five_pct_loss_to_available_end"
    else:
        reason = "five_pct_loss_to_sixty_pct_population"
    return start, end, reason


def _trajectory_metrics(per_grain: pd.DataFrame) -> tuple[float, float]:
    metrics = []
    for _, grain in per_grain.groupby("grain_id"):
        grain = grain.sort_values("time")
        if len(grain) < 3:
            continue
        time = grain["time"].to_numpy(float)
        area = grain["area"].to_numpy(float)
        metrics.append(jerkiness_metrics(time, area))
    jerk = float(np.mean([m["jerkiness_CV"] for m in metrics])) if metrics else np.nan
    burst = float(np.mean([m["burstiness"] for m in metrics])) if metrics else np.nan
    return jerk, burst


def _boundary_metrics(run_dir: Path) -> tuple[float, float, float]:
    path = run_dir / "boundary_tracks.csv"
    if not path.exists() or path.stat().st_size == 0:
        return np.nan, np.nan, np.nan
    boundaries = pd.read_csv(path)
    if boundaries.empty:
        return np.nan, np.nan, np.nan
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
    return reverse, curvature_r2, pinned


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
        renamed = radius.rename(columns={"R_A": f"R_{index}", "grain_count": f"N_{index}"})
        renamed = renamed.drop(columns="step")
        aligned = renamed if aligned is None else aligned.merge(renamed, on="time", how="inner")
    assert aligned is not None
    radius_columns = [f"R_{i}" for i in range(len(loaded))]
    count_columns = [f"N_{i}" for i in range(len(loaded))]
    radii = aligned[radius_columns].to_numpy(float).T
    mean_count = aligned[count_columns].to_numpy(float).mean(axis=1)
    start, end, window_reason = _fit_window(mean_count)
    time = aligned["time"].to_numpy(float)[start:end]
    fit_radii = radii[:, start:end]
    fit = fit_growth_law(time, fit_radii.mean(axis=0), transient_fraction=0.0)
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
        bootstrap_k.append(float(np.polyfit(time, sample_radius**fit.exponent, 1)[0]))
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
        "jerkiness_CV": float(np.nanmean([m[0] for m in trajectory_metrics])),
        "Fano": float(np.nanmean([m[1] for m in event_metrics])) if any(np.isfinite(m[1]) for m in event_metrics) else np.nan,
        "burstiness": float(np.nanmean([m[1] for m in trajectory_metrics])),
        "reverse_motion_fraction": float(np.nanmean([m[0] for m in boundary_metrics])) if any(np.isfinite(m[0]) for m in boundary_metrics) else np.nan,
        "velocity_curvature_R2": float(np.nanmean([m[1] for m in boundary_metrics])) if any(np.isfinite(m[1]) for m in boundary_metrics) else np.nan,
        "pinned_fraction": float(np.nanmean([m[2] for m in boundary_metrics])) if any(np.isfinite(m[2]) for m in boundary_metrics) else np.nan,
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
    }
    return row, diagnostics


def analyze_run(run_dir: str | Path) -> dict[str, object]:
    """Analyze one run using the same estimator as a one-member ensemble."""
    return analyze_group([Path(run_dir)], bootstrap_samples=1)[0]


def analyze_campaign(campaign_dir: str | Path, output: str | Path | None = None) -> pd.DataFrame:
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
        row, detail = analyze_group(paths)
        rows.append(row)
        diagnostics.append(detail)
    summary = pd.DataFrame(rows, columns=SUMMARY_COLUMNS)
    for regime, indices in summary.groupby("regime").groups.items():
        subset = summary.loc[indices]
        if len(subset) >= 4 and np.all(subset["K"] > 0):
            activation = fit_activation_energy(subset["temperature"].to_numpy(), subset["K"].to_numpy())
            summary.loc[indices, "Q_app"] = activation.activation_energy_ev
            summary.loc[indices, "Q_app_ci"] = 1.96 * activation.standard_error_ev
    target = Path(output) if output else campaign_dir / "mechanism_summary.csv"
    summary.to_csv(target, index=False)
    target.with_name(f"{target.stem}_diagnostics.json").write_text(json.dumps(diagnostics, indent=2))
    return summary
