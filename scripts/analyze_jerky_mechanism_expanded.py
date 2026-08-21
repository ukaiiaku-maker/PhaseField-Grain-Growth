#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from grain_growth_pf.analysis.grain_tracks import ensemble_radius, load_tracks
from grain_growth_pf.analysis.growth_law import fit_growth_law_fixed_exponent
from grain_growth_pf.analysis.jerkiness import jerkiness_metrics
from grain_growth_pf.analysis.mechanism_diagnostics import (
    activation_wait_episodes,
    arrhenius_work_summary,
    event_burst_coupling,
    factorial_effects,
    fit_pin_lifetime_models,
    kaplan_meier,
    mechanism_set,
    pin_episodes,
    reconstruct_gb_occupancy,
    release_events,
)
from grain_growth_pf.io.event_ledger import event_ledger_path, read_event_ledger


TOPOLOGY_WINDOWS = ((190, 175), (175, 160), (190, 160))


def completed_runs(root: Path) -> list[Path]:
    manifest = json.loads((root / "campaign_manifest.json").read_text())
    runs = []
    for raw in manifest["runs"]:
        run = Path(raw)
        run_manifest = json.loads((run / "manifest.json").read_text())
        if run_manifest.get("status") == "completed":
            runs.append(run)
    return runs


def load_events(run: Path) -> pd.DataFrame:
    path = event_ledger_path(run)
    if not path.exists():
        return pd.DataFrame(columns=["time", "event_type", "entity_id", "grain_ids"])
    frame = read_event_ledger(path)
    for column in ("time", "event_type", "entity_id", "grain_ids"):
        if column not in frame:
            frame[column] = pd.Series(dtype=float if column == "time" else object)
    return frame


def topology_metrics(
    tracks: pd.DataFrame,
    boundaries: pd.DataFrame,
    events: pd.DataFrame,
    high: int,
    low: int,
) -> dict[str, float]:
    radius = ensemble_radius(tracks)
    selected = radius[(radius["grain_count"] <= high) & (radius["grain_count"] >= low)]
    if len(selected) < 3:
        return {"samples": len(selected)}
    time = selected["time"].to_numpy(float)
    values = selected["R_A"].to_numpy(float)
    speed = np.abs(np.diff(values) / np.diff(time))
    tolerance = 0.05 * float(np.mean(speed)) if len(speed) else 0.0
    metrics = jerkiness_metrics(time, values, stationary_tolerance=tolerance)
    fit = fit_growth_law_fixed_exponent(time, values, 1.0, transient_fraction=0.0)
    start, end = float(time[0]), float(time[-1])
    boundary_window = boundaries[(boundaries["time"] >= start) & (boundaries["time"] <= end)]
    releases = release_events(events)
    release_times = releases.loc[(releases["time"] >= start) & (releases["time"] <= end), "time"].to_numpy(float)
    counts = np.histogram(release_times, bins=time)[0].astype(float)
    fano = float(np.var(counts) / np.mean(counts)) if len(counts) and np.mean(counts) else 0.0
    return {
        "samples": len(selected),
        "time_start": start,
        "time_end": end,
        "N_start_observed": int(selected["grain_count"].iloc[0]),
        "N_end_observed": int(selected["grain_count"].iloc[-1]),
        "R_start": float(values[0]),
        "R_end": float(values[-1]),
        "growth_rate_dR_dt": float(np.polyfit(time, values, 1)[0]),
        "K_linear_n1": fit.coefficient,
        "R2_linear_n1": fit.r_squared,
        "stationary_speed_tolerance": tolerance,
        "pinned_fraction": float(boundary_window["blocked"].mean()) if len(boundary_window) else np.nan,
        "release_count": len(release_times),
        "release_Fano": fano,
        **metrics,
    }


def pre_release_histories(
    boundaries: pd.DataFrame,
    work: pd.DataFrame,
    *,
    gb_energy: float,
    free_volume_stiffness: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    columns = [
        "event_type", "entity_id", "event_time", "episode_start", "episode_duration",
        "start_effective_DeltaG", "release_effective_DeltaG", "Delta_effective_DeltaG",
        "start_work_capillary", "release_work_capillary", "start_work_shear",
        "release_work_shear", "start_work_free_volume", "release_work_free_volume",
        "Delta_work_capillary", "Delta_work_shear", "Delta_work_free_volume",
        "history_samples",
    ]
    if boundaries.empty or work.empty:
        return pd.DataFrame(), pd.DataFrame(columns=columns)
    groups = {}
    for entity, group in boundaries.groupby("entity_id", sort=False):
        ordered = group.sort_values("time")
        groups[entity] = {
            "time": ordered["time"].to_numpy(float),
            "blocked": ordered["blocked"].to_numpy(bool),
            "curvature": ordered["curvature"].to_numpy(float),
            "shear": ordered["resolved_shear"].to_numpy(float),
            "free": ordered["free_volume_deficit"].to_numpy(float),
        }
    history_rows, event_rows = [], []
    for event in work.itertuples(index=False):
        entity = str(event.entity_id)
        if not entity.startswith("gb:") or entity not in groups:
            continue
        group = groups[entity]
        time = group["time"]
        index = int(np.searchsorted(time, float(event.time), side="right") - 1)
        if index < 0:
            continue
        gaps = np.diff(time)
        cadence = float(np.median(gaps[gaps > 0])) if np.any(gaps > 0) else 0.0
        while index >= 0 and not group["blocked"][index]:
            if float(event.time) - time[index] > 2.5 * cadence + 1e-12:
                break
            index -= 1
        if index < 0 or not group["blocked"][index]:
            continue
        start = index
        while start > 0 and group["blocked"][start - 1]:
            start -= 1
        sample_indices = np.arange(start, index + 1)
        sample_times = time[sample_indices]
        duration = max(float(event.time) - float(sample_times[0]), np.finfo(float).tiny)
        progress = (sample_times - sample_times[0]) / duration
        cap = gb_energy * group["curvature"][sample_indices] * float(event.activation_volume_normal)
        shear = group["shear"][sample_indices] * float(event.activation_volume_shear)
        free = free_volume_stiffness * group["free"][sample_indices] * float(event.activation_vacancies)
        effective = float(event.DeltaG0) - cap - shear - free
        for p, c, s, f, e in zip(progress, cap, shear, free, effective):
            history_rows.append({
                "event_type": event.event_type, "entity_id": entity,
                "event_time": float(event.time), "progress": float(np.clip(p, 0.0, 1.0)),
                "work_capillary": float(c), "work_shear": float(s),
                "work_free_volume": float(f), "effective_DeltaG": float(e),
                "source": "boundary_track",
            })
        history_rows.append({
            "event_type": event.event_type, "entity_id": entity,
            "event_time": float(event.time), "progress": 1.0,
            "work_capillary": float(event.work_capillary),
            "work_shear": float(event.work_shear),
            "work_free_volume": float(event.work_free_volume),
            "effective_DeltaG": float(event.effective_DeltaG), "source": "release_exact",
        })
        event_rows.append({
            "event_type": event.event_type, "entity_id": entity,
            "event_time": float(event.time), "episode_start": float(sample_times[0]),
            "episode_duration": duration,
            "start_effective_DeltaG": float(effective[0]),
            "release_effective_DeltaG": float(event.effective_DeltaG),
            "Delta_effective_DeltaG": float(event.effective_DeltaG - effective[0]),
            "start_work_capillary": float(cap[0]),
            "release_work_capillary": float(event.work_capillary),
            "start_work_shear": float(shear[0]),
            "release_work_shear": float(event.work_shear),
            "start_work_free_volume": float(free[0]),
            "release_work_free_volume": float(event.work_free_volume),
            "Delta_work_capillary": float(event.work_capillary - cap[0]),
            "Delta_work_shear": float(event.work_shear - shear[0]),
            "Delta_work_free_volume": float(event.work_free_volume - free[0]),
            "history_samples": len(sample_indices) + 1,
        })
    return pd.DataFrame(history_rows), pd.DataFrame(event_rows, columns=columns)


def aggregate_histories(history: pd.DataFrame, bins: int = 11) -> pd.DataFrame:
    if history.empty:
        return pd.DataFrame()
    frame = history.copy()
    frame["progress_bin"] = np.rint(frame["progress"] * (bins - 1)).astype(int)
    rows = []
    for keys, group in frame.groupby(["campaign", "regime", "event_type", "progress_bin"]):
        campaign, regime, event_type, progress_bin = keys
        row = {
            "campaign": campaign, "regime": regime, "event_type": event_type,
            "progress_bin": progress_bin, "normalized_progress": progress_bin / (bins - 1),
            "samples": len(group), "events": group["event_time"].nunique(),
        }
        for column in ("work_capillary", "work_shear", "work_free_volume", "effective_DeltaG"):
            values = group[column].to_numpy(float)
            row[f"{column}_mean"] = float(np.mean(values))
            row[f"{column}_median"] = float(np.median(values))
            row[f"{column}_p10"] = float(np.quantile(values, 0.10))
            row[f"{column}_p90"] = float(np.quantile(values, 0.90))
        rows.append(row)
    return pd.DataFrame(rows)


def tj_waiting_lower_bound(events: pd.DataFrame) -> dict[str, float]:
    durations = []
    tj = events[events["event_type"].isin({"tj_activation_hit", "tj_compatibility_release"})]
    for _, group in tj.groupby("entity_id", sort=False):
        start = None
        for row in group.sort_values("time").itertuples(index=False):
            if row.event_type == "tj_activation_hit" and start is None:
                start = float(row.time)
            elif row.event_type == "tj_compatibility_release" and start is not None:
                durations.append(float(row.time) - start)
                start = None
    values = np.asarray(durations, float)
    return {
        "tj_observed_cycles": len(values),
        "tj_post_first_hit_domain_time": float(np.sum(values)) if len(values) else 0.0,
        "tj_post_first_hit_median_wait": float(np.median(values)) if len(values) else np.nan,
        "tj_wait_scope": "lower bound from first persisted activation hit to TJ release",
    }


def minimal_model_ranking(summary: pd.DataFrame) -> pd.DataFrame:
    factorial = summary[summary["regime"].map(lambda value: value == "B0" or set(value) <= set("GTSC"))].copy()
    factorial = factorial[factorial["regime"] != "QIU"]
    metrics = [
        "K_linear_n1", "stationary_fraction", "motion_top_5pct", "Fano",
        "median_pin_duration", "top_5pct_release_excess",
    ]
    factorial[["Fano", "median_pin_duration", "top_5pct_release_excess"]] = factorial[
        ["Fano", "median_pin_duration", "top_5pct_release_excess"]
    ].fillna(0.0)
    target = factorial.loc[factorial["regime"] == "GTSC", metrics].iloc[0].to_numpy(float)
    values = factorial[metrics].to_numpy(float)
    scale = np.std(values, axis=0, ddof=0)
    scale[scale <= np.finfo(float).tiny] = 1.0
    factorial["GTSC_standardized_RMS_discrepancy"] = np.sqrt(np.mean(((values - target) / scale) ** 2, axis=1))
    factorial["mechanism_count"] = factorial["regime"].map(lambda value: len(mechanism_set(value)))
    pareto = []
    for row in factorial.itertuples():
        dominated = np.any(
            (factorial["mechanism_count"] <= row.mechanism_count)
            & (factorial["GTSC_standardized_RMS_discrepancy"] <= row.GTSC_standardized_RMS_discrepancy)
            & ((factorial["mechanism_count"] < row.mechanism_count)
               | (factorial["GTSC_standardized_RMS_discrepancy"] < row.GTSC_standardized_RMS_discrepancy))
        )
        pareto.append(not bool(dominated))
    factorial["Pareto_optimal"] = pareto
    return factorial[["regime", "mechanism_count", "GTSC_standardized_RMS_discrepancy", "Pareto_optimal", *metrics]].sort_values(
        ["mechanism_count", "GTSC_standardized_RMS_discrepancy"]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("factorial_root")
    parser.add_argument("--discriminating-root")
    parser.add_argument("--factorial-summary", required=True)
    parser.add_argument("--output-dir", default="results/production_summaries/jerky_mechanism_expanded")
    parser.add_argument("--shuffle-samples", type=int, default=200)
    args = parser.parse_args()

    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    campaigns = [("factorial", Path(args.factorial_root))]
    if args.discriminating_root:
        campaigns.append(("discriminating", Path(args.discriminating_root)))

    causal_rows, topology_rows, survival_fit_rows, survival_curve_rows = [], [], [], []
    multiplier_rows, history_frames, history_event_frames, occupancy_rows = [], [], [], []
    for campaign_name, root in campaigns:
        for run_index, run in enumerate(completed_runs(root)):
            manifest = json.loads((run / "manifest.json").read_text())
            config = manifest["config"]
            regime = config["regime"]
            print(f"{campaign_name}: {regime}", flush=True)
            tracks = load_tracks(run / "grain_tracks.csv")
            boundaries = pd.read_csv(run / "boundary_tracks.csv")
            events = load_events(run)
            work_path = run / "activation_work.csv"
            work = pd.read_csv(work_path) if work_path.exists() else pd.DataFrame()

            causal = event_burst_coupling(
                tracks, events, shuffle_samples=args.shuffle_samples,
                seed=20260821 + run_index + (1000 if campaign_name == "discriminating" else 0),
            )
            causal_rows.append({"campaign": campaign_name, "regime": regime, "run": str(run), **causal})

            episodes = pin_episodes(boundaries)
            fitted = fit_pin_lifetime_models(episodes)
            if not fitted.empty:
                fitted.insert(0, "scope", "gb_blocked_track")
                fitted.insert(0, "regime", regime)
                fitted.insert(0, "campaign", campaign_name)
                survival_fit_rows.append(fitted)
            curve = kaplan_meier(episodes)
            if not curve.empty:
                curve.insert(0, "scope", "gb_blocked_track")
                curve.insert(0, "regime", regime)
                curve.insert(0, "campaign", campaign_name)
                survival_curve_rows.append(curve)

            tj_episodes = activation_wait_episodes(
                events,
                hit_type="tj_activation_hit",
                release_type="tj_compatibility_release",
                end_time=float(tracks["time"].max()),
            )
            tj_fitted = fit_pin_lifetime_models(tj_episodes)
            if not tj_fitted.empty:
                tj_fitted.insert(0, "scope", "tj_post_first_hit_lower_bound")
                tj_fitted.insert(0, "regime", regime)
                tj_fitted.insert(0, "campaign", campaign_name)
                survival_fit_rows.append(tj_fitted)
            tj_curve = kaplan_meier(tj_episodes)
            if not tj_curve.empty:
                tj_curve.insert(0, "scope", "tj_post_first_hit_lower_bound")
                tj_curve.insert(0, "regime", regime)
                tj_curve.insert(0, "campaign", campaign_name)
                survival_curve_rows.append(tj_curve)

            active = set(config["active_modules"])
            occupancy_rows.append({
                "campaign": campaign_name, "regime": regime,
                **reconstruct_gb_occupancy(
                    episodes, events,
                    has_g=bool(active.intersection({"gb_compatibility", "gb_area_point_defect_pinning", "gb_pinning"})),
                    has_c=bool(active.intersection({"free_volume", "serial_climb", "nucleation_limited", "exchange_limited", "transport_limited"})),
                ),
                **tj_waiting_lower_bound(events),
            })

            if not work.empty:
                multipliers = arrhenius_work_summary(work, config["pf"]["temperature"])
                multipliers.insert(0, "regime", regime)
                multipliers.insert(0, "campaign", campaign_name)
                multiplier_rows.append(multipliers)
                history, event_changes = pre_release_histories(
                    boundaries, work, gb_energy=float(config["pf"]["gb_energy"]),
                    free_volume_stiffness=float(config["parameters"].get("free_volume_stiffness", 0.0)),
                )
                if not history.empty:
                    history.insert(0, "regime", regime)
                    history.insert(0, "campaign", campaign_name)
                    history_frames.append(history)
                if not event_changes.empty:
                    event_changes.insert(0, "regime", regime)
                    event_changes.insert(0, "campaign", campaign_name)
                    history_event_frames.append(event_changes)

            if campaign_name == "factorial":
                for high, low in TOPOLOGY_WINDOWS:
                    topology_rows.append({
                        "campaign": campaign_name, "regime": regime,
                        "window": f"N{high}_to_{low}", "N_high": high, "N_low": low,
                        **topology_metrics(tracks, boundaries, events, high, low),
                    })

    causal_table = pd.DataFrame(causal_rows)
    causal_table.to_csv(output / "event_burst_causality.csv", index=False)
    pd.DataFrame(topology_rows).to_csv(output / "common_topology_windows.csv", index=False)
    pd.concat(survival_fit_rows, ignore_index=True).to_csv(output / "pin_lifetime_model_fits.csv", index=False)
    pd.concat(survival_curve_rows, ignore_index=True).to_csv(output / "pin_lifetime_survival.csv", index=False)
    pd.concat(multiplier_rows, ignore_index=True).to_csv(output / "arrhenius_work_multipliers.csv", index=False)
    pd.DataFrame(occupancy_rows).to_csv(output / "mechanism_occupancy_reconstruction.csv", index=False)

    histories = pd.concat(history_frames, ignore_index=True) if history_frames else pd.DataFrame()
    aggregate_histories(histories).to_csv(output / "pre_release_barrier_trajectories.csv", index=False)
    changes = pd.concat(history_event_frames, ignore_index=True) if history_event_frames else pd.DataFrame()
    changes.to_csv(output / "pre_release_event_changes.csv", index=False)

    summary = pd.read_csv(args.factorial_summary)
    summary = summary.merge(
        causal_table[causal_table["campaign"] == "factorial"].drop(columns=["campaign", "run"]),
        on="regime", how="left",
    )
    factorial_metrics = [
        "jerkiness_CV", "stationary_fraction", "motion_top_1pct", "motion_top_5pct", "motion_top_10pct",
        "Fano", "median_pin_duration", "K_linear_n1", "large_burst_risk_ratio_excess",
        "top_1pct_release_excess", "top_5pct_release_excess", "top_10pct_release_excess",
        "event_growth_xcorr_excess",
    ]
    design = summary[summary["regime"].map(lambda value: value == "B0" or (set(value) <= set("GTSC") and value != "QIU"))].copy()
    design[["Fano", "median_pin_duration"]] = design[["Fano", "median_pin_duration"]].fillna(0.0)
    effects = factorial_effects(design, factorial_metrics)
    effects.to_csv(output / "factorial_effects.csv", index=False)
    ranking = minimal_model_ranking(summary)
    ranking.to_csv(output / "minimal_model_ranking.csv", index=False)

    report = [
        "# Expanded jerky-mechanism analysis", "",
        f"Factorial campaign: `{args.factorial_root}` (16 G/T/S/C subsets plus QIU = 17 cases).",
        f"Discriminating campaign: `{args.discriminating_root}`." if args.discriminating_root else "",
        "", "No new phase-field simulations were launched.", "",
        "## Interpretation guardrails", "",
        "- Event-burst coupling uses grain-linked releases and one recorded-frame lookback. Its shuffled control preserves each grain's number of event-associated intervals.",
        "- Factorial effects are descriptive one-realization contrasts; they have no seed-level uncertainty estimate.",
        "- G/C occupancy is event-classified from the generic GB blocked flag plus completion ledgers. Episodes with no persisted completion remain unresolved. Exact G/T/C reason-specific occupancy was not persisted; TJ waiting is only a lower bound from first saved hit to release. Shear is a hazard modifier, not a gate.",
        "- Pre-release histories use saved boundary frames and append the exact release work. Episodes shorter than the boundary-output cadence cannot be reconstructed.",
        "- Lifetime fits include right censoring in the likelihood. AIC compares exponential (one parameter), gamma, and Weibull (two parameters).", "",
        "## Outputs", "",
    ]
    for path in sorted(output.glob("*.csv")):
        report.append(f"- `{path.name}`")
    (output / "README.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(output / "README.md")


if __name__ == "__main__":
    main()
