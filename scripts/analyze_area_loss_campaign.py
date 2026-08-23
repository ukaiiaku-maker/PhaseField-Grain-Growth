#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp

from grain_growth_pf.io.event_ledger import event_ledger_path, read_event_ledger


COMPLETIONS = {"gb_sink_completion", "tj_sink_completion"}
WINDOWS = (1, 2, 4, 8)


def _runs(root: Path) -> list[Path]:
    for name in ("campaign_manifest.json", "video_manifest.json"):
        path = root / name
        if path.exists():
            manifest = json.loads(path.read_text())
            return [
                Path(item["path"] if isinstance(item, dict) else item)
                for item in manifest["runs"]
            ]
    raise ValueError(f"no campaign/video manifest under {root}")


def _events(run: Path) -> pd.DataFrame:
    path = event_ledger_path(run)
    if not path.exists():
        return pd.DataFrame()
    return read_event_ledger(path)


def _position(value: Any) -> np.ndarray | None:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    if isinstance(value, str):
        try:
            value = json.loads(value)
        except json.JSONDecodeError:
            return None
    result = np.asarray(value, dtype=float)
    return result if result.shape == (2,) and np.all(np.isfinite(result)) else None


def _growth_series(tracks: pd.DataFrame) -> pd.DataFrame:
    grouped = tracks.groupby(["step", "time"], as_index=False).agg(
        mean_radius=("radius", "mean"), grain_count=("grain_id", "nunique")
    ).sort_values("step")
    if len(grouped) > 1:
        grouped["radius_rate"] = np.gradient(
            grouped["mean_radius"].to_numpy(float), grouped["time"].to_numpy(float)
        )
    else:
        grouped["radius_rate"] = 0.0
    return grouped


def _stage_audit(regime: str, events: pd.DataFrame) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    stage_rows: list[dict[str, Any]] = []
    cycle_rows: list[dict[str, Any]] = []
    if events.empty or "event_type" not in events:
        return stage_rows, cycle_rows
    for path in ("gb_disconnection_sink", "tj_sink"):
        selected = events[events["event_type"].isin({
            f"{path}_nucleation", f"{path}_exchange", f"{path}_transport"
        })].copy()
        for stage in ("nucleation", "exchange", "transport"):
            data = selected[selected["event_type"] == f"{path}_{stage}"]
            residence = pd.to_numeric(data.get("stage_residence_time"), errors="coerce").dropna()
            rate = pd.to_numeric(data.get("instantaneous_rate"), errors="coerce")
            expected = (1.0 / rate[rate > 0.0]).replace([np.inf, -np.inf], np.nan).dropna()
            stage_rows.append({
                "regime": regime, "sink_path": path, "stage": stage,
                "events": len(data),
                "observed_mean_residence": residence.mean() if len(residence) else np.nan,
                "observed_median_residence": residence.median() if len(residence) else np.nan,
                "rate_implied_mean_residence": expected.mean() if len(expected) else np.nan,
                "median_rate": rate.median() if len(rate) else np.nan,
            })
        complete_type = "gb_sink_completion" if path.startswith("gb_") else "tj_sink_completion"
        complete = events[events["event_type"] == complete_type].sort_values(["entity_id", "time"])
        cycle_times: list[float] = []
        cycle_rates: list[tuple[float, float, float]] = []
        for _, entity_rows in selected.sort_values(["entity_id", "time"]).groupby("entity_id"):
            accumulated = 0.0
            seen: list[str] = []
            rates: list[float] = []
            for _, event in entity_rows.iterrows():
                stage = str(event["event_type"]).removeprefix(path + "_")
                residence = float(event.get("stage_residence_time", np.nan))
                rate = float(event.get("instantaneous_rate", np.nan))
                if np.isfinite(residence):
                    accumulated += residence
                    seen.append(stage)
                    rates.append(rate)
                if stage == "transport":
                    if seen == ["nucleation", "exchange", "transport"]:
                        cycle_times.append(accumulated)
                        cycle_rates.append(tuple(rates))
                    accumulated = 0.0
                    seen = []
                    rates = []
        implied_means = [
            row["rate_implied_mean_residence"] for row in stage_rows[-3:]
            if np.isfinite(row["rate_implied_mean_residence"])
        ]
        predicted: list[float] = []
        if cycle_rates:
            rng = np.random.default_rng(20260822)
            draws_per_cycle = max(16, min(256, 100_000 // len(cycle_rates)))
            for rates in cycle_rates:
                if all(np.isfinite(rate) and rate > 0.0 for rate in rates):
                    sample = sum(rng.exponential(1.0 / rate, draws_per_cycle) for rate in rates)
                    predicted.extend(sample.tolist())
        ks = ks_2samp(cycle_times, predicted) if cycle_times and predicted else None
        cycle_rows.append({
            "regime": regime, "sink_path": path, "completed_cycles": len(complete),
            "observed_cycle_intervals": len(cycle_times),
            "observed_mean_cycle_time": np.mean(cycle_times) if cycle_times else np.nan,
            "observed_median_cycle_time": np.median(cycle_times) if cycle_times else np.nan,
            "observed_p10_cycle_time": np.quantile(cycle_times, 0.1) if cycle_times else np.nan,
            "observed_p90_cycle_time": np.quantile(cycle_times, 0.9) if cycle_times else np.nan,
            "serial_stage_mean_prediction": sum(implied_means) if len(implied_means) == 3 else np.nan,
            "hypoexponential_predicted_mean": np.mean(predicted) if predicted else np.nan,
            "hypoexponential_predicted_median": np.median(predicted) if predicted else np.nan,
            "hypoexponential_predicted_p10": np.quantile(predicted, 0.1) if predicted else np.nan,
            "hypoexponential_predicted_p90": np.quantile(predicted, 0.9) if predicted else np.nan,
            "hypoexponential_ks_statistic": ks.statistic if ks is not None else np.nan,
            "hypoexponential_ks_pvalue": ks.pvalue if ks is not None else np.nan,
            "hypoexponential_model": "conditional_stage_rate_mixture",
        })
    return stage_rows, cycle_rows


def _local_responses(regime: str, boundary: pd.DataFrame, events: pd.DataFrame) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if boundary.empty or events.empty:
        return rows
    releases = events[events["event_type"].isin(COMPLETIONS)]
    grouped = {key: frame.sort_values("step").reset_index(drop=True) for key, frame in boundary.groupby("entity_id")}
    for _, event in releases.iterrows():
        candidates: list[pd.DataFrame] = []
        if event["entity_id"] in grouped:
            candidates.append(grouped[event["entity_id"]])
        elif event["event_type"] == "tj_sink_completion":
            grains = {int(item) for item in str(event.get("grain_ids", "")).split(";") if item}
            for frame in grouped.values():
                if {int(frame.iloc[0]["grain_i"]), int(frame.iloc[0]["grain_j"])} <= grains:
                    candidates.append(frame)
        for window in WINDOWS:
            signed: list[float] = []
            absolute: list[float] = []
            displacement: list[float] = []
            for frame in candidates:
                index = int(np.argmin(np.abs(frame["step"].to_numpy(float) - float(event["step"]))))
                if index < window or index + window >= len(frame):
                    continue
                before = float(frame.iloc[index - window]["normal_velocity"])
                after = float(frame.iloc[index + window]["normal_velocity"])
                signed.append(after - before)
                absolute.append(abs(after) - abs(before))
                local = frame.iloc[index - window:index + window + 1]
                displacement.append(float(np.trapezoid(
                    local["normal_velocity"].to_numpy(float), local["time"].to_numpy(float)
                )))
            if signed:
                rows.append({
                    "regime": regime, "event_type": event["event_type"],
                    "event_id": event.get("event_id"), "window_frames": window,
                    "domains": len(signed), "signed_velocity_response": np.mean(signed),
                    "absolute_velocity_response": np.mean(absolute),
                    "signed_normal_displacement": np.mean(displacement),
                })
    return rows


def _causal_nulls(regime: str, growth: pd.DataFrame, events: pd.DataFrame) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if len(growth) < 10 or events.empty or "event_type" not in events:
        return rows
    steps = growth["step"].to_numpy(int)
    response = np.abs(growth["radius_rate"].to_numpy(float))
    indicator = np.zeros(len(steps), dtype=bool)
    for event_step in events.loc[events["event_type"].isin(COMPLETIONS), "step"].astype(int):
        insertion = int(np.searchsorted(steps, event_step))
        candidates = [index for index in (insertion - 1, insertion) if 0 <= index < len(steps)]
        if candidates:
            indicator[min(candidates, key=lambda index: abs(int(steps[index]) - event_step))] = True
    rng = np.random.default_rng(20260822)
    for window in WINDOWS:
        def score(flags: np.ndarray) -> float:
            indices = np.flatnonzero(flags)
            values = [response[i + 1:i + 1 + window].mean() for i in indices if i + window < len(response)]
            return float(np.mean(values)) if values else np.nan

        actual = score(indicator)
        circular = []
        block = []
        block_size = max(8, window)
        blocks = [indicator[i:i + block_size] for i in range(0, len(indicator), block_size)]
        for _ in range(200):
            shift = int(rng.integers(1, len(indicator)))
            circular.append(score(np.roll(indicator, shift)))
            order = rng.permutation(len(blocks))
            shuffled = np.concatenate([blocks[i] for i in order])[:len(indicator)]
            block.append(score(shuffled))
        circular_mean = float(np.mean([value for value in circular if np.isfinite(value)])) if any(np.isfinite(circular)) else np.nan
        block_mean = float(np.mean([value for value in block if np.isfinite(value)])) if any(np.isfinite(block)) else np.nan
        rows.append({
            "regime": regime, "window_frames": window,
            "actual_abs_motion": actual,
            "circular_null_mean": circular_mean,
            "circular_excess": actual - circular_mean,
            "block_null_mean": block_mean,
            "block_excess": actual - block_mean,
            "release_frames": int(indicator.sum()), "shuffles": 200,
        })
    return rows


def _resistance_windows(regime: str, growth: pd.DataFrame) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if len(growth) < 12:
        return rows
    data = growth[np.isfinite(growth["radius_rate"]) & (growth["radius_rate"] > 0.0)].copy()
    if len(data) < 8:
        return rows
    data["window"] = pd.qcut(np.arange(len(data)), q=min(5, len(data)), labels=False, duplicates="drop")
    for window, frame in data.groupby("window"):
        if len(frame) < 3:
            continue
        x = frame["mean_radius"].to_numpy(float)
        y = 1.0 / frame["radius_rate"].to_numpy(float)
        slope, intercept = np.polyfit(x, y, 1)
        prediction = slope * x + intercept
        ss_total = float(np.sum((y - y.mean()) ** 2))
        rows.append({
            "regime": regime, "topology_window": int(window),
            "grain_count_start": int(frame["grain_count"].iloc[0]),
            "grain_count_end": int(frame["grain_count"].iloc[-1]),
            "points": len(frame), "a": slope, "b": intercept,
            "r2": 1.0 - float(np.sum((y - prediction) ** 2)) / ss_total if ss_total else np.nan,
        })
    return rows


def _release_selection(
    regime: str,
    boundary: pd.DataFrame,
    events: pd.DataFrame,
    inventory: pd.DataFrame,
    free_volume_stiffness: float,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if boundary.empty or events.empty:
        return rows
    rng = np.random.default_rng(5101)
    grouped = {key: frame.sort_values("step") for key, frame in boundary.groupby("entity_id")}
    for _, event in events[events["event_type"] == "gb_sink_completion"].iterrows():
        frame = grouped.get(event["entity_id"])
        if frame is None or frame.empty:
            continue
        release = frame.iloc[int(np.argmin(np.abs(frame["step"].to_numpy(float) - float(event["step"]))))]
        pinned = frame[frame["blocked"] > 0]
        control = pinned.iloc[int(rng.integers(len(pinned)))] if len(pinned) else frame.iloc[int(rng.integers(len(frame)))]
        quota = abs(float(event.get("signed_defect_quota", 0.0) or 0.0))
        for sample, state in (("release", release), ("matched_pinned", control)):
            if sample == "release":
                chemical_work = float(event.get("work_chemical", np.nan))
                work_total = float(event.get("work_total", np.nan))
            elif not inventory.empty:
                index = int(np.argmin(np.abs(
                    inventory["step"].to_numpy(float) - float(state["step"])
                )))
                stored_signed = float(inventory.iloc[index]["N_stored_signed"])
                chemical_work = float(free_volume_stiffness * stored_signed * quota)
                work_total = chemical_work
            else:
                chemical_work = work_total = np.nan
            rows.append({
                "regime": regime, "event_id": event.get("event_id"), "sample": sample,
                "pVn": float(state["curvature"] * state["normal_velocity"]),
                "tau_Vtau": float(state["resolved_shear"] * abs(state["normal_velocity"])),
                "DeltaMu_Nv": chemical_work,
                "W_total": work_total,
            })
    return rows


def _occupancy(regime: str, state: pd.DataFrame) -> list[dict[str, Any]]:
    rows = []
    if state.empty:
        return rows
    for entity_type, frame in state.groupby("entity_type"):
        weights = frame.groupby("step")["time"].first().sort_index().diff().shift(-1)
        default = float(weights.dropna().median()) if weights.notna().any() else 1.0
        weights = weights.fillna(default).to_dict()
        w = frame["step"].map(weights).to_numpy(float)
        denominator = float(w.sum())
        rows.append({
            "regime": regime, "entity_type": entity_type,
            "entity_time": denominator,
            "G_fraction": float(np.sum(w * frame["G_pending"]) / denominator),
            "T_fraction": float(np.sum(w * frame["T_pending"]) / denominator),
            "C_fraction": float(np.sum(w * frame["C_pending"]) / denominator),
            "blocked_fraction": float(np.sum(w * frame["blocked"]) / denominator),
            "G_and_C_fraction": float(np.sum(w * (frame["G_pending"] * frame["C_pending"])) / denominator),
        })
    return rows


def _grain_set(value: Any) -> frozenset[int]:
    return frozenset(int(item) for item in str(value).split(";") if item and item != "nan")


def _gb_arclength_coordinates(boundary: pd.DataFrame) -> dict[str, tuple[tuple[int, int], float]]:
    """Estimate persistent-domain centers along each GB from saved domain lengths."""
    if boundary.empty:
        return {}
    lengths = boundary.groupby(["grain_i", "grain_j", "entity_id"], as_index=False)["length"].median()
    coordinates: dict[str, tuple[tuple[int, int], float]] = {}
    for (grain_i, grain_j), frame in lengths.groupby(["grain_i", "grain_j"]):
        indexed: list[tuple[int, str, float]] = []
        for row in frame.itertuples():
            try:
                index = int(str(row.entity_id).rsplit(":", 1)[1])
            except (IndexError, ValueError):
                continue
            indexed.append((index, str(row.entity_id), float(row.length)))
        cursor = 0.0
        for _, entity_id, length in sorted(indexed):
            coordinates[entity_id] = ((int(grain_i), int(grain_j)), cursor + 0.5 * length)
            cursor += length
    return coordinates


def _shares_tj(first: frozenset[int], second: frozenset[int]) -> bool:
    if len(first) == 2 and len(second) == 2:
        return len(first & second) == 1 and len(first | second) == 3
    if len(first) == 2 and len(second) == 3:
        return first <= second
    if len(first) == 3 and len(second) == 2:
        return second <= first
    if len(first) == 3 and len(second) == 3:
        return len(first & second) >= 2
    return False


def _spatial_correlations(
    regime: str, events: pd.DataFrame, boundary: pd.DataFrame, shape: tuple[int, int]
) -> list[dict[str, Any]]:
    if events.empty or "event_type" not in events:
        return []
    releases = events[events["event_type"].isin(COMPLETIONS)].sort_values("time")
    records = []
    arclength = _gb_arclength_coordinates(boundary)
    parsed = [(_position(row.position), row) for row in releases.itertuples()]
    max_lag = 1.0
    max_pairs = 200_000
    for i, (p0, first) in enumerate(parsed):
        if p0 is None:
            continue
        for p1, second in parsed[i + 1:]:
            lag = float(second.time - first.time)
            if lag > max_lag:
                break
            if p1 is None:
                continue
            delta = np.abs(p1 - p0)
            box = np.asarray(shape, dtype=float)
            delta = np.minimum(delta, box - delta)
            distance = float(np.linalg.norm(delta))
            grains0 = _grain_set(first.grain_ids)
            grains1 = _grain_set(second.grain_ids)
            first_arc = arclength.get(str(first.entity_id))
            second_arc = arclength.get(str(second.entity_id))
            along_gb = np.nan
            if first_arc is not None and second_arc is not None and first_arc[0] == second_arc[0]:
                along_gb = abs(first_arc[1] - second_arc[1])
            records.append({
                "regime": regime, "time_lag": lag, "distance": distance,
                "same_entity": int(first.entity_id == second.entity_id),
                "gb_arclength_separation": along_gb,
                "same_gb": int(first_arc is not None and second_arc is not None and first_arc[0] == second_arc[0]),
                "shared_tj_connectivity": int(_shares_tj(grains0, grains1)),
                "first_sink": first.event_type, "second_sink": second.event_type,
            })
            if len(records) >= max_pairs:
                return records
    return records


def _spatial_summary(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    frame = pd.DataFrame(records)
    if frame.empty:
        return []
    rows: list[dict[str, Any]] = []
    time_bins = pd.cut(frame["time_lag"], [-np.inf, 0.04, 0.16, 0.32, 1.0], include_lowest=True)
    for metric, bins in (
        ("distance", [-np.inf, 4.0, 8.0, 16.0, 32.0, 64.0, np.inf]),
        ("gb_arclength_separation", [-np.inf, 6.0, 12.0, 24.0, 48.0, np.inf]),
    ):
        valid = frame[np.isfinite(frame[metric])].copy()
        if valid.empty:
            continue
        valid["metric_bin"] = pd.cut(valid[metric], bins, include_lowest=True).astype(str)
        valid["time_bin"] = time_bins.loc[valid.index].astype(str)
        for (regime, metric_bin, time_bin), group in valid.groupby(
            ["regime", "metric_bin", "time_bin"], observed=True
        ):
            rows.append({
                "regime": regime, "metric": metric, "metric_bin": metric_bin,
                "time_bin": time_bin, "pairs": len(group),
                "same_entity_fraction": group["same_entity"].mean(),
                "same_gb_fraction": group["same_gb"].mean(),
                "shared_tj_connectivity_fraction": group["shared_tj_connectivity"].mean(),
            })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    root = Path(args.campaign)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)

    tables: dict[str, list[dict[str, Any]]] = {
        name: [] for name in (
            "run_summary", "stage_residence", "cycle_time_audit", "local_event_response",
            "causal_nulls", "kinetic_resistance_windows", "release_state_selection",
            "direct_occupancy", "sink_partition", "spatial_release_correlations",
            "spatial_release_summary",
        )
    }
    for run in _runs(root):
        manifest = json.loads((run / "manifest.json").read_text())
        config = manifest["config"]
        regime = str(config["regime"])
        tracks = pd.read_csv(run / "grain_tracks.csv")
        boundary = pd.read_csv(run / "boundary_tracks.csv")
        events = _events(run)
        state_path = run / "mechanism_state.csv"
        state = pd.read_csv(state_path) if state_path.exists() else pd.DataFrame()
        inventory_path = run / "defect_inventory.csv"
        inventory = pd.read_csv(inventory_path) if inventory_path.exists() else pd.DataFrame()
        growth = _growth_series(tracks)
        start = manifest.get("restart_provenance", {})
        tables["run_summary"].append({
            "regime": regime, "path": str(run), "source_sha": manifest["git_sha"],
            "starting_step": int(start.get("checkpoint_step", 0)),
            "ending_step": int(manifest["steps_completed"]),
            "starting_time": float(start.get("checkpoint_time", 0.0)),
            "ending_time": float(json.loads((run / "checkpoint.json").read_text())["time"]),
            "starting_grains": int(start.get("starting_grains", manifest.get("grains_after_equilibration", growth["grain_count"].iloc[0]))),
            "ending_grains": int(manifest["final_grains"]),
            "status": manifest["status"],
        })
        stage, cycles = _stage_audit(regime, events)
        tables["stage_residence"].extend(stage)
        tables["cycle_time_audit"].extend(cycles)
        tables["local_event_response"].extend(_local_responses(regime, boundary, events))
        tables["causal_nulls"].extend(_causal_nulls(regime, growth, events))
        tables["kinetic_resistance_windows"].extend(_resistance_windows(regime, growth))
        tables["release_state_selection"].extend(_release_selection(
            regime, boundary, events, inventory,
            float(config["parameters"].get("free_volume_stiffness", 0.05)),
        ))
        tables["direct_occupancy"].extend(_occupancy(regime, state))
        if not inventory.empty:
            final = inventory.iloc[-1]
            tables["sink_partition"].append({
                "regime": regime, "N_required": final["N_required"],
                "N_accommodated_GB": final["N_accommodated_GB"],
                "N_accommodated_TJ": final["N_accommodated_TJ"],
                "f_GB": final["GB_accommodation_fraction"],
                "f_TJ": final["TJ_accommodation_fraction"],
                "stored": final["N_active_deficit"] + final["N_retired"],
                "conservation_residual": final["conservation_residual"],
                "max_abs_conservation_residual": final["max_abs_conservation_residual"],
            })
        tables["spatial_release_correlations"].extend(
            _spatial_correlations(regime, events, boundary, tuple(config["pf"]["shape"]))
        )

    tables["spatial_release_summary"] = _spatial_summary(
        tables["spatial_release_correlations"]
    )

    for name, rows in tables.items():
        pd.DataFrame(rows).to_csv(output / f"{name}.csv", index=False)
    (output / "analysis_manifest.json").write_text(json.dumps({
        "campaign": str(root), "outputs": {
            name: str(output / f"{name}.csv") for name in tables
        }, "response_windows_frames": WINDOWS,
        "causal_nulls": ["circular_time_shift", "block_shuffle"],
    }, indent=2) + "\n")
    print(output)


if __name__ == "__main__":
    main()
