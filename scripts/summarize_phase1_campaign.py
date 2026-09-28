#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np
import pandas as pd


REGIME_ORDER = (
    "B0", "GTSC_GBTJ_Ks025", "GTSC_GBTJ_Ks030_LONG", "GTC_GB",
    "GSC_GBTJ_Ks025", "G", "T", "GT", "QIU",
)

CLASSIFICATIONS = {
    "B0": "apparently converged long-time growth regime",
    "G": "apparently converged long-time growth regime",
    "T": "apparently converged long-time growth regime",
    "GT": "apparently converged long-time growth regime",
    "GTC_GB": "apparently converged long-time growth regime",
    "GSC_GBTJ_Ks025": "progressive slowing over the available dynamic range",
    "GTSC_GBTJ_Ks025": "progressive slowing over the available dynamic range",
    "GTSC_GBTJ_Ks030_LONG": "progressive slowing over the available dynamic range",
    "QIU": "transient still evolving toward a regime",
}

INTERPRETATIONS = {
    "B0": "finite exponent near 1.8 cumulatively; late local fit is finite-size noisy",
    "G": "finite exponent near 1.5 cumulatively; G mainly reduces the coefficient",
    "T": "finite exponent near 1.5 cumulatively; no persistent high-n branch",
    "GT": "finite exponent near 1.5 cumulatively; no persistent high-n branch",
    "GTC_GB": "finite exponent near 1.6 cumulatively with stable late K2",
    "GSC_GBTJ_Ks025": "n rises through 7.3, 11, 23, and 45 while K2 collapses",
    "GTSC_GBTJ_Ks025": "n rises through 7.2, 9.9, and 30 while K2 collapses",
    "GTSC_GBTJ_Ks030_LONG": "n reaches the search bound (>=50) and terminal G is flat",
    "QIU": "late avalanche from N=494 to 99 invalidates a single limiting-law claim",
}

NEXT_ACTION = {
    "B0": "none; terminal N=100 criterion reached",
    "G": "none; terminal N=100 criterion reached",
    "T": "none; terminal N=100 criterion reached",
    "GT": "none; terminal N=100 criterion reached",
    "GTC_GB": "none; terminal N=100 criterion reached",
    "GSC_GBTJ_Ks025": "checkpoint-continue to 150000 steps, then reapply the kinetic gate",
    "GTSC_GBTJ_Ks025": "checkpoint-continue to 150000 steps, then reapply the kinetic gate",
    "GTSC_GBTJ_Ks030_LONG": "retain as progressive-slowing result; extend only to test persistence, not to estimate n",
    "QIU": "diagnose the step-10000 avalanche and morphology before any continuation",
}


def _sha256(path: Path) -> str:
    import hashlib

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _crossing_time(trajectory: pd.DataFrame, target: float) -> float:
    x = trajectory.G_population_over_G0.to_numpy(float)
    time = trajectory.time.to_numpy(float)
    crossing = np.flatnonzero(x >= target)
    if not len(crossing):
        return float("nan")
    index = int(crossing[0])
    if index == 0 or x[index] == x[index - 1]:
        return float(time[index])
    fraction = (target - x[index - 1]) / (x[index] - x[index - 1])
    return float(time[index - 1] + fraction * (time[index] - time[index - 1]))


def _common_metrics(trajectory: pd.DataFrame, lo: float = 1.15, hi: float = 1.40) -> dict:
    ta, tb = _crossing_time(trajectory, lo), _crossing_time(trajectory, hi)
    g0 = float(trajectory.G_population.iloc[0])
    ga, gb = lo * g0, hi * g0
    dt = tb - ta
    if not np.isfinite(dt) or dt <= 0:
        return {"available": False}
    return {
        "available": True,
        "x_start": lo,
        "x_end": hi,
        "time_start": ta,
        "time_end": tb,
        "elapsed": dt,
        "G_rate": (gb - ga) / dt,
        "K2_secant": (gb**2 - ga**2) / dt,
        "K3_secant": (gb**3 - ga**3) / dt,
    }


def _direct_occupancy(run: Path, regime: str) -> list[dict]:
    path = run / "mechanism_state.csv"
    if not path.is_file():
        return []
    sums: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for chunk in pd.read_csv(path, chunksize=250_000):
        for entity_type, group in chunk.groupby("entity_type"):
            bucket = sums[str(entity_type)]
            bucket["entity_samples"] += len(group)
            for source, target in (
                ("G_pending", "G_fraction"), ("T_pending", "T_fraction"),
                ("C_pending", "C_fraction"), ("blocked", "blocked_fraction"),
            ):
                bucket[target] += pd.to_numeric(group[source], errors="coerce").fillna(0).sum()
            bucket["G_and_C_fraction"] += (
                (pd.to_numeric(group.G_pending, errors="coerce").fillna(0) > 0)
                & (pd.to_numeric(group.C_pending, errors="coerce").fillna(0) > 0)
            ).sum()
    rows = []
    for entity_type, values in sorted(sums.items()):
        count = values["entity_samples"]
        rows.append(
            {
                "regime": regime,
                "entity_type": entity_type,
                "entity_samples": int(count),
                **{
                    name: values[name] / count
                    for name in (
                        "G_fraction", "T_fraction", "C_fraction",
                        "blocked_fraction", "G_and_C_fraction",
                    )
                },
            }
        )
    return rows


def _stage_rows(run: Path, regime: str) -> list[dict]:
    event_root = run / "events.parquet"
    if not list(event_root.glob("*.parquet")):
        return []
    frame = pd.read_parquet(
        event_root, columns=["event_type", "sink_path", "stage_residence_time"]
    )
    frame = frame[frame.stage_residence_time.notna()]
    rows = []
    for (sink, event), group in frame.groupby(["sink_path", "event_type"]):
        values = pd.to_numeric(group.stage_residence_time, errors="coerce").dropna()
        if not len(values) or not any(stage in str(event) for stage in ("nucleation", "exchange", "transport")):
            continue
        rows.append(
            {
                "regime": regime,
                "sink_path": sink,
                "stage": str(event).rsplit("_", 1)[-1],
                "events": len(values),
                "mean_residence": float(values.mean()),
                "median_residence": float(values.median()),
                "p90_residence": float(values.quantile(0.9)),
                "p99_residence": float(values.quantile(0.99)),
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    args = parser.parse_args()
    root = Path(args.campaign)
    analysis = root / "analysis"
    campaign = json.loads((root / "campaign_manifest.json").read_text(encoding="utf-8"))
    integration = json.loads((root / "integration_manifest.json").read_text(encoding="utf-8"))
    trajectories = pd.read_csv(analysis / "grain_size_trajectory.csv")
    fits = pd.read_csv(analysis / "growth_law_windows.csv")
    collapse = pd.read_csv(analysis / "distribution_collapse.csv")
    summary = pd.read_csv(analysis / "long_time_run_summary.csv").set_index("regime")
    run_by_regime = {}
    for raw in campaign["runs"]:
        run = Path(raw)
        manifest = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
        run_by_regime[manifest["config"]["regime"]] = run

    classifications, common, sink_rows, occupancy_rows, stage_rows, movie_rows = [], [], [], [], [], []
    for regime in REGIME_ORDER:
        run = run_by_regime[regime]
        trajectory = trajectories[trajectories.regime == regime].sort_values("step")
        local_fits = fits[
            (fits.regime == regime) & (fits.size_measure == "G_population")
            & (fits.available.fillna(False)) & fits.window.str.startswith("x_")
        ].sort_values("window_end_x")
        cumulative = fits[
            (fits.regime == regime) & (fits.size_measure == "G_population")
            & (fits.available.fillna(False)) & fits.window.str.startswith("cumulative")
        ].sort_values("window_end_x")
        grain_changes = trajectory[trajectory.grain_count.diff().fillna(0) != 0]
        last_change_step = int(grain_changes.step.iloc[-1]) if len(grain_changes) else 0
        latest_fit = local_fits.iloc[-1] if len(local_fits) else None
        cumulative_fit = cumulative.iloc[-1] if len(cumulative) else None
        latest_collapse = collapse[
            (collapse.regime == regime) & (collapse.metric == "equivalent_diameter")
        ].sort_values("actual_x")
        classifications.append(
            {
                "regime": regime,
                "classification": CLASSIFICATIONS[regime],
                "interpretation": INTERPRETATIONS[regime],
                "next_action": NEXT_ACTION[regime],
                "end_step": int(summary.loc[regime, "steps"]),
                "end_N": int(summary.loc[regime, "end_N"]),
                "end_G_over_G0": float(summary.loc[regime, "G_over_G0"]),
                "target_N100_reached": bool(summary.loc[regime, "target_N100_reached"]),
                "last_grain_count_change_step": last_change_step,
                "terminal_steps_since_grain_count_change": int(summary.loc[regime, "steps"]) - last_change_step,
                "latest_window": latest_fit.window if latest_fit is not None else "",
                "latest_n_best": float(latest_fit.n_best) if latest_fit is not None else np.nan,
                "latest_n_profile_width_5pct": float(latest_fit.n_profile_width_5pct) if latest_fit is not None else np.nan,
                "latest_K2": float(latest_fit.K2) if latest_fit is not None else np.nan,
                "latest_K3": float(latest_fit.K3) if latest_fit is not None else np.nan,
                "cumulative_n_best": float(cumulative_fit.n_best) if cumulative_fit is not None else np.nan,
                "cumulative_K2": float(cumulative_fit.K2) if cumulative_fit is not None else np.nan,
                "latest_size_KS": float(latest_collapse.ks_distance.iloc[-1]) if len(latest_collapse) else np.nan,
                "terminal_mean_compactness": float(summary.loc[regime, "terminal_mean_compactness"]),
                "max_abs_conservation_residual": float(summary.loc[regime, "max_abs_conservation_residual"]),
            }
        )
        common.append({"regime": regime, **_common_metrics(trajectory)})
        occupancy_rows.extend(_direct_occupancy(run, regime))
        stage_rows.extend(_stage_rows(run, regime))
        inventory_path = run / "defect_inventory.csv"
        if inventory_path.is_file():
            inventory = pd.read_csv(inventory_path)
            terminal = inventory.iloc[-1]
            sink_rows.append(
                {
                    "regime": regime,
                    "N_required": terminal.N_required,
                    "N_accommodated_GB": terminal.N_accommodated_GB,
                    "N_accommodated_TJ": terminal.N_accommodated_TJ,
                    "N_stored": terminal.N_active_deficit + terminal.N_retired,
                    "GB_accommodation_fraction": terminal.GB_accommodation_fraction,
                    "TJ_accommodation_fraction": terminal.TJ_accommodation_fraction,
                    "max_abs_conservation_residual": inventory.max_abs_conservation_residual.max(),
                }
            )
        movie = run / "microstructure.mp4"
        movie_rows.append(
            {
                "regime": regime,
                "movie_path": str(movie),
                "movie_exists": movie.is_file(),
                "movie_size_bytes": movie.stat().st_size if movie.is_file() else 0,
                "movie_sha256": _sha256(movie) if movie.is_file() else "",
                "source_frame_count": len(list((run / "frames").glob("*.npz"))),
                "first_step": int(trajectory.step.iloc[0]),
                "last_step": int(trajectory.step.iloc[-1]),
            }
        )

    common_frame = pd.DataFrame(common).set_index("regime")
    contrasts = {
        "G minus B0": {"G": 1, "B0": -1},
        "T minus B0": {"T": 1, "B0": -1},
        "G:T interaction": {"GT": 1, "G": -1, "T": -1, "B0": 1},
        "C_GB added to GT": {"GTC_GB": 1, "GT": -1},
        "T added to GSC_GBTJ_Ks025": {"GTSC_GBTJ_Ks025": 1, "GSC_GBTJ_Ks025": -1},
        "full GBTJ plus S versus GTC_GB": {"GTSC_GBTJ_Ks025": 1, "GTC_GB": -1},
        "Ks 0.30 minus 0.25 full": {"GTSC_GBTJ_Ks030_LONG": 1, "GTSC_GBTJ_Ks025": -1},
    }
    factorial_rows = []
    for name, weights in contrasts.items():
        for metric in ("G_rate", "K2_secant", "K3_secant", "elapsed"):
            factorial_rows.append(
                {
                    "common_window": "x=1.15 to 1.40",
                    "contrast": name,
                    "metric": metric,
                    "contrast_value": sum(weight * common_frame.loc[regime, metric] for regime, weight in weights.items()),
                    "terms": ";".join(f"{weight:+g}*{regime}" for regime, weight in weights.items()),
                }
            )

    outputs = {
        "phase1_regime_classification.csv": pd.DataFrame(classifications),
        "phase1_common_progress_metrics.csv": pd.DataFrame(common),
        "phase1_revised_factorial_comparisons.csv": pd.DataFrame(factorial_rows),
        "phase1_sink_partition.csv": pd.DataFrame(sink_rows),
        "phase1_direct_occupancy.csv": pd.DataFrame(occupancy_rows),
        "phase1_stage_resolved_kinetics.csv": pd.DataFrame(stage_rows),
        "phase1_movie_index.csv": pd.DataFrame(movie_rows),
    }
    for name, frame in outputs.items():
        frame.to_csv(analysis / name, index=False)
    payload = {
        "schema_version": 1,
        "campaign_root": str(root),
        "campaign_status": campaign["status"],
        "all_nine_integrated": integration["all_nine_integrated"],
        "scientific_source_commit": campaign["git_sha"],
        "initial_state_npz_sha256": integration["initial_state_npz_sha256"],
        "analysis_tables": {name: len(frame) for name, frame in outputs.items()},
        "movies_complete": all(row["movie_exists"] for row in movie_rows),
    }
    causal = analysis / "causal_nulls_vs_grain_size.csv"
    if causal.is_file():
        payload["analysis_tables"][causal.name] = len(pd.read_csv(causal))
        payload["causal_null_sha256"] = _sha256(causal)
    payload["contact_sheets"] = {
        name: {
            "exists": (analysis / name).is_file(),
            "sha256": _sha256(analysis / name) if (analysis / name).is_file() else "",
        }
        for name in ("terminal_contact_sheet.png", "progress_matched_contact_sheet.png")
    }
    payload["artifact_inventory"] = [
        {
            "path": path.name,
            "size_bytes": path.stat().st_size,
            "sha256": _sha256(path),
        }
        for path in sorted(analysis.iterdir())
        if path.is_file() and path.name != "phase1_analysis_manifest.json"
    ]
    (analysis / "phase1_analysis_manifest.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(analysis)


if __name__ == "__main__":
    main()
