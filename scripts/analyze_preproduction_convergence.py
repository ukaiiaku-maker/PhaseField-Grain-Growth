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


def _runs(campaign: Path) -> list[Path]:
    manifest = json.loads((campaign / "campaign_manifest.json").read_text())
    result = []
    for raw in manifest.get("runs", []):
        run = Path(raw)
        mpath = run / "manifest.json"
        if not mpath.exists():
            continue
        data = json.loads(mpath.read_text())
        if data.get("status") == "completed":
            result.append(run)
    if not result:
        raise ValueError(f"no completed runs in {campaign}")
    return result


def _mechanism(regime: str) -> str:
    for name in ("B0", "T2", "C5"):
        if regime.startswith(name):
            return name
    return regime.split("_")[0]


def _boundary_metrics(path: Path) -> tuple[float, float, float]:
    frame = pd.read_csv(path / "boundary_tracks.csv")
    if frame.empty:
        return np.nan, np.nan, np.nan
    k = pd.to_numeric(frame["curvature"], errors="coerce").to_numpy(float)
    v = pd.to_numeric(frame["normal_velocity"], errors="coerce").to_numpy(float)
    valid = np.isfinite(k) & np.isfinite(v) & (np.abs(k) > 1e-12) & (np.abs(v) > 1e-12)
    if np.any(valid):
        kt = np.quantile(np.abs(k[valid]), 0.75)
        vt = np.quantile(np.abs(v[valid]), 0.75)
        active = valid & (np.abs(k) >= kt) & (np.abs(v) >= vt)
        reverse = float(np.mean(k[active] * v[active] < 0)) if np.any(active) else np.nan
    else:
        reverse = np.nan
    if np.count_nonzero(valid) > 2 and np.std(k[valid]) and np.std(v[valid]):
        corr = np.corrcoef(k[valid], v[valid])[0, 1]
        r2 = float(corr * corr) if np.isfinite(corr) else np.nan
    else:
        r2 = np.nan
    pinned = float(np.mean(pd.to_numeric(frame["blocked"], errors="coerce").fillna(0.0)))
    return reverse, r2, pinned


def _jerkiness(per_grain: pd.DataFrame) -> float:
    values = []
    for _, grain in per_grain.groupby("grain_id"):
        grain = grain.sort_values("time")
        if len(grain) < 3:
            continue
        metrics = jerkiness_metrics(
            grain["time"].to_numpy(float), grain["area"].to_numpy(float)
        )
        value = metrics.get("jerkiness_CV", np.nan)
        if np.isfinite(value):
            values.append(float(value))
    return float(np.mean(values)) if values else np.nan


def _time_to_fraction(radius: pd.DataFrame, initial_n: float, fraction: float) -> float:
    hit = radius[radius["grain_count"] <= fraction * initial_n]
    return float(hit.iloc[0]["time"]) if not hit.empty else np.nan


def _run_metrics(run: Path, scenario: str) -> dict[str, object]:
    manifest = json.loads((run / "manifest.json").read_text())
    config = manifest["config"]
    regime = str(config["regime"])
    per_grain = load_tracks(run / "grain_tracks.csv")
    radius = ensemble_radius(per_grain)
    initial_n = float(radius.iloc[0]["grain_count"])
    selection = (
        (radius["grain_count"] <= 0.95 * initial_n)
        & (radius["grain_count"] >= 0.50 * initial_n)
    )
    fit_frame = radius.loc[selection]
    if len(fit_frame) < 8:
        fit_frame = radius.iloc[max(1, int(0.05 * len(radius))):]
    fit = fit_growth_law_fixed_exponent(
        fit_frame["time"].to_numpy(float),
        fit_frame["R_A"].to_numpy(float),
        2.0,
        transient_fraction=0.0,
    )

    q_frame = per_grain[(per_grain["area"] > 0) & np.isfinite(per_grain["perimeter"])].copy()
    q_frame["q"] = q_frame["perimeter"] ** 2 / (4.0 * np.pi * q_frame["area"])
    q_by_time = q_frame.groupby("time", as_index=False)["q"].mean()
    q_window = q_by_time[
        (q_by_time["time"] >= float(fit_frame["time"].min()))
        & (q_by_time["time"] <= float(fit_frame["time"].max()))
    ]
    q_mean = float(q_window["q"].mean()) if not q_window.empty else np.nan

    reverse, curvature_r2, pinned = _boundary_metrics(run)
    return {
        "scenario": scenario,
        "regime": regime,
        "mechanism": _mechanism(regime),
        "shape_y": int(config["pf"]["shape"][0]),
        "shape_x": int(config["pf"]["shape"][1]),
        "dx": float(config["pf"]["grid_spacing"]),
        "dt": float(config["pf"]["time_step"]),
        "physical_size_y": float(config["pf"]["shape"][0]) * float(config["pf"]["grid_spacing"]),
        "initial_grains": initial_n,
        "final_grains": float(radius.iloc[-1]["grain_count"]),
        "K2": float(fit.coefficient),
        "K2_r2": float(fit.r_squared),
        "q_mean": q_mean,
        "jerkiness_CV": _jerkiness(per_grain),
        "reverse_motion_fraction": reverse,
        "velocity_curvature_R2": curvature_r2,
        "pinned_fraction": pinned,
        "t_N75pct": _time_to_fraction(radius, initial_n, 0.75),
        "t_N625pct": _time_to_fraction(radius, initial_n, 0.625),
        "t_N50pct": _time_to_fraction(radius, initial_n, 0.50),
        "run": str(run),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True)
    parser.add_argument("--dt-half", required=True)
    parser.add_argument("--grid-fine", required=True)
    parser.add_argument("--size-large", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    campaigns = {
        "reference": Path(args.reference),
        "dt_half": Path(args.dt_half),
        "grid_fine": Path(args.grid_fine),
        "size_large": Path(args.size_large),
    }
    rows = []
    for scenario, campaign in campaigns.items():
        rows.extend(_run_metrics(run, scenario) for run in _runs(campaign))
    table = pd.DataFrame(rows).sort_values(["scenario", "mechanism"]).reset_index(drop=True)

    for scenario in table["scenario"].unique():
        idx = table["scenario"] == scenario
        subset = table.loc[idx]
        b0 = subset[subset["mechanism"] == "B0"].iloc[0]
        table.loc[idx, "K2_over_B0"] = subset["K2"].to_numpy(float) / float(b0["K2"])
        table.loc[idx, "q_over_B0"] = subset["q_mean"].to_numpy(float) / float(b0["q_mean"])
        for field in ("t_N75pct", "t_N625pct", "t_N50pct"):
            denom = float(b0[field])
            table.loc[idx, field + "_over_B0"] = (
                subset[field].to_numpy(float) / denom if np.isfinite(denom) and denom > 0 else np.nan
            )

    ref = table[table["scenario"] == "reference"].set_index("mechanism")
    for i, row in table.iterrows():
        mechanism = row["mechanism"]
        if mechanism not in ref.index:
            continue
        r = ref.loc[mechanism]
        table.loc[i, "K2_change_vs_reference"] = row["K2"] / r["K2"] - 1.0
        table.loc[i, "K2_over_B0_change_vs_reference"] = row["K2_over_B0"] / r["K2_over_B0"] - 1.0
        table.loc[i, "q_over_B0_change_vs_reference"] = row["q_over_B0"] / r["q_over_B0"] - 1.0

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, index=False)

    report = output.with_suffix(".md")
    lines = [
        "# Preproduction numerical-convergence check",
        "",
        "Single-realization artifact screen. This is not a statistical replication campaign.",
        "Primary invariant is each mechanism's K2/K2(B0) and morphology relative to B0.",
        "",
    ]
    for scenario in ("reference", "dt_half", "grid_fine", "size_large"):
        lines.append(f"## {scenario}")
        subset = table[table["scenario"] == scenario]
        for _, row in subset.iterrows():
            lines.append(
                f"- {row['mechanism']}: K2={row['K2']:.6g}, "
                f"K2/B0={row['K2_over_B0']:.4f}, q/B0={row['q_over_B0']:.4f}, "
                f"jerkiness={row['jerkiness_CV']:.4f}, reverse={row['reverse_motion_fraction']:.4f}"
            )
        if scenario != "reference":
            active = subset[subset["mechanism"].isin(["T2", "C5"])]
            if len(active):
                max_k = float(np.nanmax(np.abs(active["K2_over_B0_change_vs_reference"])))
                max_q = float(np.nanmax(np.abs(active["q_over_B0_change_vs_reference"])))
                lines.append(f"- max |change in K2/B0| vs reference: {max_k:.2%}")
                lines.append(f"- max |change in q/B0| vs reference: {max_q:.2%}")
        lines.append("")
    lines += [
        "## Interpretation gate",
        "",
        "Treat approximately <=15% change in normalized kinetic resistance and <=5% change in normalized compactness as a practical preproduction screen, not a fitted uncertainty interval. Large systematic changes require diagnosis before model freeze.",
        "",
    ]
    report.write_text("\n".join(lines))
    print(table.to_string(index=False))
    print(report)


if __name__ == "__main__":
    main()
