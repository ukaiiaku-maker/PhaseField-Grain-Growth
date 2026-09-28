#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


TARGET_COUNTS = (150, 125, 100)
FIT_START_COUNT = 190
FIT_END_COUNT = 100
BOOTSTRAP_SAMPLES = 20000
RNG_SEED = 20260819


def _run_dirs(campaign: Path) -> list[Path]:
    manifest = json.loads((campaign / "campaign_manifest.json").read_text())
    raw = manifest.get("runs", [])
    result = []
    for item in raw:
        path = Path(item)
        if not path.is_absolute() and not path.exists():
            path = campaign.parent.parent / path
        result.append(path)
    return result


def _observables(run_dir: Path) -> tuple[dict, pd.DataFrame]:
    manifest = json.loads((run_dir / "manifest.json").read_text())
    tracks = pd.read_csv(run_dir / "grain_tracks.csv")
    tracks = tracks[(tracks["area"] > 0) & (tracks["perimeter"] > 0)].copy()
    tracks["compactness_q"] = (
        tracks["perimeter"].to_numpy(float) ** 2
        / (4.0 * np.pi * tracks["area"].to_numpy(float))
    )
    grouped = tracks.groupby(["time", "step"], as_index=False).agg(
        grain_count=("grain_id", "nunique"),
        mean_area=("area", "mean"),
        q_mean=("compactness_q", "mean"),
        q_median=("compactness_q", "median"),
        q_p90=("compactness_q", lambda x: float(np.quantile(x, 0.9))),
    )
    grouped["R_A"] = np.sqrt(grouped["mean_area"] / np.pi)
    return manifest, grouped.sort_values(["time", "step"]).reset_index(drop=True)


def _first_at_or_below(frame: pd.DataFrame, count: int) -> pd.Series | None:
    rows = frame[frame["grain_count"] <= count]
    return None if rows.empty else rows.iloc[0]


def _fixed_n2_slope(frame: pd.DataFrame) -> tuple[float, int, float, float]:
    start_rows = frame[frame["grain_count"] <= FIT_START_COUNT]
    if start_rows.empty:
        start_index = 0
    else:
        start_index = int(start_rows.index[0])
    end_rows = frame[(frame.index >= start_index) & (frame["grain_count"] <= FIT_END_COUNT)]
    end_index = int(end_rows.index[0]) if not end_rows.empty else int(frame.index[-1])
    window = frame.loc[start_index:end_index].copy()
    if len(window) < 8:
        return np.nan, len(window), np.nan, np.nan
    t = window["time"].to_numpy(float)
    r2 = window["R_A"].to_numpy(float) ** 2
    dt = t - t[0]
    dr2 = r2 - r2[0]
    denominator = float(dt @ dt)
    if denominator <= 0:
        return np.nan, len(window), float(window["grain_count"].iloc[0]), float(window["grain_count"].iloc[-1])
    slope = float(dt @ dr2 / denominator)
    return slope, len(window), float(window["grain_count"].iloc[0]), float(window["grain_count"].iloc[-1])


def _bootstrap_mean_ci(values: np.ndarray) -> tuple[float, float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan, np.nan, np.nan
    if len(values) == 1:
        value = float(values[0])
        return value, value, value
    rng = np.random.default_rng(RNG_SEED)
    samples = rng.choice(values, size=(BOOTSTRAP_SAMPLES, len(values)), replace=True).mean(axis=1)
    return float(values.mean()), float(np.quantile(samples, 0.025)), float(np.quantile(samples, 0.975))


def _per_run(campaign: Path) -> pd.DataFrame:
    records = []
    for run_dir in _run_dirs(campaign):
        manifest, obs = _observables(run_dir)
        config = manifest["config"]
        slope, samples, n_start, n_end = _fixed_n2_slope(obs)
        record = {
            "regime": config["regime"],
            "seed": int(config["seed"]),
            "K2_190_to_100": slope,
            "fit_samples": samples,
            "fit_start_grains": n_start,
            "fit_end_grains": n_end,
            "final_grains": int(obs["grain_count"].iloc[-1]),
            "final_time": float(obs["time"].iloc[-1]),
        }
        t0 = float(obs["time"].iloc[0])
        for count in TARGET_COUNTS:
            row = _first_at_or_below(obs, count)
            record[f"time_to_N{count}"] = np.nan if row is None else float(row["time"] - t0)
            record[f"q_mean_N{count}"] = np.nan if row is None else float(row["q_mean"])
            record[f"q_p90_N{count}"] = np.nan if row is None else float(row["q_p90"])
        records.append(record)
    return pd.DataFrame(records).sort_values(["regime", "seed"]).reset_index(drop=True)


def _paired_to_b0(per_run: pd.DataFrame) -> pd.DataFrame:
    baseline = per_run[per_run["regime"] == "B0_CTL"].set_index("seed")
    rows = []
    for _, row in per_run[per_run["regime"] != "B0_CTL"].iterrows():
        seed = int(row["seed"])
        if seed not in baseline.index:
            continue
        base = baseline.loc[seed]
        result = {"regime": row["regime"], "seed": seed}
        for column in ["K2_190_to_100", *[f"time_to_N{x}" for x in TARGET_COUNTS], *[f"q_mean_N{x}" for x in TARGET_COUNTS]]:
            numerator = float(row[column])
            denominator = float(base[column])
            result[f"{column}_over_B0"] = (
                numerator / denominator
                if np.isfinite(numerator) and np.isfinite(denominator) and denominator != 0
                else np.nan
            )
        rows.append(result)
    return pd.DataFrame(rows)


def _grouped(per_run: pd.DataFrame, paired: pd.DataFrame, standard_summary: Path | None) -> pd.DataFrame:
    rows = []
    standard = None
    if standard_summary is not None and standard_summary.exists():
        standard = pd.read_csv(standard_summary).set_index("regime")
    for regime, group in per_run.groupby("regime", sort=True):
        row = {
            "regime": regime,
            "realizations": len(group),
            "K2_mean": float(group["K2_190_to_100"].mean()),
            "K2_sd": float(group["K2_190_to_100"].std(ddof=1)) if len(group) > 1 else np.nan,
            "final_grains_mean": float(group["final_grains"].mean()),
        }
        if regime != "B0_CTL" and not paired.empty:
            p = paired[paired["regime"] == regime]
            for column in ["K2_190_to_100", *[f"time_to_N{x}" for x in TARGET_COUNTS], *[f"q_mean_N{x}" for x in TARGET_COUNTS]]:
                metric = f"{column}_over_B0"
                mean, low, high = _bootstrap_mean_ci(p[metric].to_numpy(float))
                row[f"{metric}_mean"] = mean
                row[f"{metric}_ci_low"] = low
                row[f"{metric}_ci_high"] = high
        for count in TARGET_COUNTS:
            row[f"time_to_N{count}_mean"] = float(group[f"time_to_N{count}"].mean())
            row[f"q_mean_N{count}_mean"] = float(group[f"q_mean_N{count}"].mean())
        if standard is not None and regime in standard.index:
            s = standard.loc[regime]
            for column in (
                "jerkiness_CV", "Fano", "burstiness", "reverse_motion_fraction",
                "velocity_curvature_R2", "pinned_fraction", "number_of_events",
            ):
                if column in s:
                    row[column] = s[column]
        rows.append(row)
    return pd.DataFrame(rows)


def _sensitivity(per_run: pd.DataFrame) -> dict[str, float]:
    line = per_run[per_run["regime"] == "C5_GATE_LINE"].set_index("seed")
    diffuse = per_run[per_run["regime"] == "C5_GATE_DIFFUSE"].set_index("seed")
    seeds = sorted(set(line.index) & set(diffuse.index))
    if not seeds:
        return {}
    k_ratio = np.asarray([diffuse.loc[s, "K2_190_to_100"] / line.loc[s, "K2_190_to_100"] for s in seeds], float)
    q_ratio = np.asarray([diffuse.loc[s, "q_mean_N100"] / line.loc[s, "q_mean_N100"] for s in seeds], float)
    km, kl, kh = _bootstrap_mean_ci(k_ratio)
    qm, ql, qh = _bootstrap_mean_ci(q_ratio)
    return {
        "C5_diffuse_over_line_K2_mean": km,
        "C5_diffuse_over_line_K2_ci_low": kl,
        "C5_diffuse_over_line_K2_ci_high": kh,
        "C5_diffuse_over_line_q_N100_mean": qm,
        "C5_diffuse_over_line_q_N100_ci_low": ql,
        "C5_diffuse_over_line_q_N100_ci_high": qh,
    }


def _report(grouped: pd.DataFrame, sensitivity: dict[str, float], output: Path) -> None:
    activated = ["G2_GATE_LINE", "T2_GATE_LINE", "S2_GATE_LINE", "C5_GATE_LINE"]
    lines = [
        "# Five-seed migration-closure validation",
        "",
        f"Primary kinetic fit: fixed n=2 over the common population window N={FIT_START_COUNT} to N={FIT_END_COUNT}.",
        "Free-n fits from the standard campaign summary are secondary diagnostics.",
        "",
        "## Primary decision metrics",
        "",
        "| regime | K2/B0 mean [95% CI] | t(N=100)/B0 | q(N=100)/B0 | jerkiness CV | reverse fraction | v-kappa R2 |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for regime in activated:
        rows = grouped[grouped["regime"] == regime]
        if rows.empty:
            continue
        r = rows.iloc[0]
        k = r.get("K2_190_to_100_over_B0_mean", np.nan)
        kl = r.get("K2_190_to_100_over_B0_ci_low", np.nan)
        kh = r.get("K2_190_to_100_over_B0_ci_high", np.nan)
        t = r.get("time_to_N100_over_B0_mean", np.nan)
        q = r.get("q_mean_N100_over_B0_mean", np.nan)
        lines.append(
            f"| {regime} | {k:.3f} [{kl:.3f}, {kh:.3f}] | {t:.3f} | {q:.3f} | "
            f"{r.get('jerkiness_CV', np.nan):.3f} | {r.get('reverse_motion_fraction', np.nan):.3f} | "
            f"{r.get('velocity_curvature_R2', np.nan):.3f} |"
        )
    lines += ["", "## Predeclared checks", ""]
    for regime in activated:
        rows = grouped[grouped["regime"] == regime]
        if rows.empty:
            continue
        r = rows.iloc[0]
        kh = r.get("K2_190_to_100_over_B0_ci_high", np.nan)
        qh = r.get("q_mean_N100_over_B0_ci_high", np.nan)
        kinetic = "PASS" if np.isfinite(kh) and kh < 1.05 else "REVIEW"
        morphology = "PASS" if np.isfinite(qh) and qh < 1.05 else "REVIEW"
        lines.append(f"- {regime}: kinetic resistance {kinetic}; matched-population morphology {morphology}.")
    if sensitivity:
        lines += [
            "",
            "## C5 line/diffuse sensitivity",
            "",
            f"K2 diffuse/line = {sensitivity['C5_diffuse_over_line_K2_mean']:.3f} "
            f"[{sensitivity['C5_diffuse_over_line_K2_ci_low']:.3f}, {sensitivity['C5_diffuse_over_line_K2_ci_high']:.3f}]",
            f"q(N=100) diffuse/line = {sensitivity['C5_diffuse_over_line_q_N100_mean']:.3f} "
            f"[{sensitivity['C5_diffuse_over_line_q_N100_ci_low']:.3f}, {sensitivity['C5_diffuse_over_line_q_N100_ci_high']:.3f}]",
        ]
    lines += [
        "",
        "PASS/REVIEW labels are screening aids, not hypothesis-test p-values. The production decision should use the full paired-seed distributions and physical trends.",
    ]
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze the five-seed gate-only migration-closure validation campaign.")
    parser.add_argument("campaign")
    parser.add_argument("--standard-summary")
    parser.add_argument("--output-prefix", default="results/production_summaries/migration_closure_validation_5seed")
    args = parser.parse_args()

    campaign = Path(args.campaign)
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    standard = Path(args.standard_summary) if args.standard_summary else None

    per_run = _per_run(campaign)
    paired = _paired_to_b0(per_run)
    grouped = _grouped(per_run, paired, standard)
    sensitivity = _sensitivity(per_run)

    per_run.to_csv(prefix.with_name(prefix.name + "_per_run.csv"), index=False)
    paired.to_csv(prefix.with_name(prefix.name + "_paired_to_B0.csv"), index=False)
    grouped.to_csv(prefix.with_name(prefix.name + "_grouped.csv"), index=False)
    _report(grouped, sensitivity, prefix.with_name(prefix.name + "_decision.md"))

    print(grouped.to_string(index=False))
    print(prefix.with_name(prefix.name + "_decision.md"))


if __name__ == "__main__":
    main()
