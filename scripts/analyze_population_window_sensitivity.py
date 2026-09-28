#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from grain_growth_pf.analysis.grain_tracks import ensemble_radius, load_tracks
from grain_growth_pf.analysis.growth_law import fit_growth_law, fit_growth_law_fixed_exponent

BANDS = ((0.95, 0.75), (0.75, 0.625), (0.625, 0.50))


def _mechanism(regime: str) -> str:
    for name in ("B0", "G2", "T2", "S2", "C5"):
        if regime.startswith(name):
            return name
    return regime


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    campaign = Path(args.campaign)
    manifest = json.loads((campaign / "campaign_manifest.json").read_text())
    rows = []
    for raw in manifest.get("runs", []):
        run = Path(raw)
        mpath = run / "manifest.json"
        if not mpath.exists():
            continue
        rmanifest = json.loads(mpath.read_text())
        if rmanifest.get("status") != "completed":
            continue
        regime = str(rmanifest["config"]["regime"])
        mechanism = _mechanism(regime)
        # Use the corrected gate-only cases and B0. Ignore legacy hybrid controls.
        if "HYBRID" in regime or "DIFFUSE" in regime:
            continue
        tracks = load_tracks(run / "grain_tracks.csv")
        radius = ensemble_radius(tracks)
        initial_n = float(radius.iloc[0]["grain_count"])
        for upper, lower in BANDS:
            selected = radius[
                (radius["grain_count"] <= upper * initial_n)
                & (radius["grain_count"] >= lower * initial_n)
            ]
            if len(selected) < 8:
                continue
            time = selected["time"].to_numpy(float)
            r = selected["R_A"].to_numpy(float)
            fixed = fit_growth_law_fixed_exponent(time, r, 2.0, transient_fraction=0.0)
            free = fit_growth_law(time, r, transient_fraction=0.0)
            rows.append({
                "regime": regime,
                "mechanism": mechanism,
                "upper_fraction": upper,
                "lower_fraction": lower,
                "N_start_approx": upper * initial_n,
                "N_end_approx": lower * initial_n,
                "samples": len(selected),
                "K2": float(fixed.coefficient),
                "K2_r2": float(fixed.r_squared),
                "free_n": float(free.exponent),
                "free_n_r2": float(free.r_squared),
                "free_n_at_bound": bool(free.exponent <= 1.01 or free.exponent >= 5.99),
                "run": str(run),
            })
    table = pd.DataFrame(rows)
    if table.empty:
        raise SystemExit("no analyzable corrected closure runs found")

    for (upper, lower), indices in table.groupby(["upper_fraction", "lower_fraction"]).groups.items():
        subset = table.loc[indices]
        b0 = subset[subset["mechanism"] == "B0"]
        if b0.empty:
            continue
        k0 = float(b0.iloc[0]["K2"])
        table.loc[indices, "K2_over_B0"] = subset["K2"].to_numpy(float) / k0

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, index=False)

    report = output.with_suffix(".md")
    lines = [
        "# Population-window sensitivity",
        "",
        "Single-trajectory analysis of the already-generated corrected closure cases.",
        "The fixed n=2 coefficient is primary; free-n fits are diagnostic only.",
        "",
    ]
    for mechanism in sorted(table["mechanism"].unique()):
        lines.append(f"## {mechanism}")
        sub = table[table["mechanism"] == mechanism].sort_values("upper_fraction", ascending=False)
        for _, row in sub.iterrows():
            lines.append(
                f"- N/N0={row['upper_fraction']:.3f}->{row['lower_fraction']:.3f}: "
                f"K2={row['K2']:.6g}, K2/B0={row.get('K2_over_B0', np.nan):.4f}, "
                f"free n={row['free_n']:.3f}{' [BOUND]' if row['free_n_at_bound'] else ''}"
            )
        if len(sub) > 1:
            k = sub["K2"].to_numpy(float)
            spread = (np.nanmax(k) - np.nanmin(k)) / max(np.nanmean(k), np.finfo(float).tiny)
            lines.append(f"- relative K2 spread across population bands: {spread:.2%}")
        lines.append("")
    report.write_text("\n".join(lines))
    print(table.to_string(index=False))
    print(report)


if __name__ == "__main__":
    main()
