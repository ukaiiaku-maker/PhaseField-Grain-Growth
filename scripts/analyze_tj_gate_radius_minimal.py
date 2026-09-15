#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from grain_growth_pf.analysis.grain_tracks import ensemble_radius, load_tracks
from grain_growth_pf.analysis.growth_law import fit_growth_law_fixed_exponent


def _campaign_runs(root: Path) -> list[Path]:
    manifest = json.loads((root / "campaign_manifest.json").read_text())
    return [Path(path) for path in manifest.get("runs", [])]


def _video_runs(root: Path) -> list[Path]:
    manifest = json.loads((root / "video_manifest.json").read_text())
    runs = manifest.get("runs", [])
    result = []
    for item in runs:
        result.append(Path(item["path"] if isinstance(item, dict) else item))
    return result


def _find(runs: list[Path], prefix: str) -> Path:
    matches = [path for path in runs if path.name.startswith(prefix)]
    if len(matches) != 1:
        raise ValueError(f"expected exactly one {prefix!r} run, found {matches}")
    return matches[0]


def _time_to_count(radius: pd.DataFrame, target: int) -> float:
    hit = radius[radius["grain_count"] <= target]
    return float(hit.iloc[0]["time"]) if not hit.empty else np.nan


def _metrics(path: Path, label: str) -> dict[str, float | str]:
    manifest = json.loads((path / "manifest.json").read_text())
    config = manifest["config"]
    tracks = load_tracks(path / "grain_tracks.csv")
    radius = ensemble_radius(tracks)
    window = radius[(radius["grain_count"] <= 190) & (radius["grain_count"] >= 100)]
    if len(window) < 8:
        raise ValueError(f"insufficient N=190..100 samples for {path}: {len(window)}")
    fit = fit_growth_law_fixed_exponent(
        window["time"].to_numpy(float), window["R_A"].to_numpy(float),
        2.0, transient_fraction=0.0,
    )
    boundary_path = path / "boundary_tracks.csv"
    pinned = np.nan
    if boundary_path.exists():
        boundary = pd.read_csv(boundary_path)
        selected = boundary[
            (boundary["time"] >= float(window["time"].iloc[0]))
            & (boundary["time"] <= float(window["time"].iloc[-1]))
        ]
        if not selected.empty:
            pinned = float(selected["blocked"].mean())
    dx = float(config["pf"]["grid_spacing"])
    radius_pixels = int(config.get("parameters", {}).get("tj_correlation_radius", 2))
    width = (2 * radius_pixels + 1) * dx
    area = (2 * radius_pixels + 1) ** 2 * dx**2
    return {
        "label": label,
        "regime": str(config["regime"]),
        "dx": dx,
        "tj_radius_pixels": radius_pixels,
        "tj_gate_width_approx": width,
        "tj_gate_area_approx": area,
        "K2": float(fit.coefficient),
        "K2_r2": float(fit.r_squared),
        "pinned_fraction_190_100": pinned,
        "t_N150": _time_to_count(radius, 150),
        "t_N125": _time_to_count(radius, 125),
        "t_N100": _time_to_count(radius, 100),
        "run": str(path),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare legacy and corrected fine-grid T2 TJ pinning footprints."
    )
    parser.add_argument("--reference-campaign", required=True)
    parser.add_argument("--legacy-fine-campaign", required=True)
    parser.add_argument("--corrected-video-root", required=True)
    parser.add_argument(
        "--output", default="results/production_summaries/tj_gate_radius_minimal.csv"
    )
    args = parser.parse_args()

    reference = _campaign_runs(Path(args.reference_campaign))
    legacy_fine = _campaign_runs(Path(args.legacy_fine_campaign))
    corrected = _video_runs(Path(args.corrected_video_root))

    rows = [
        _metrics(_find(reference, "B0_REF-"), "B0_reference_dx1"),
        _metrics(_find(reference, "T2_REF-"), "T2_reference_dx1_r2"),
        _metrics(_find(legacy_fine, "B0_GRID_FINE-"), "B0_fine_dx075"),
        _metrics(_find(legacy_fine, "T2_GRID_FINE-"), "T2_fine_legacy_r2"),
        _metrics(_find(corrected, "T2_GRID_FINE_TJ_R3-"), "T2_fine_corrected_r3"),
    ]
    frame = pd.DataFrame(rows)

    ref_b0 = float(frame.loc[frame["label"] == "B0_reference_dx1", "K2"].iloc[0])
    ref_t2 = float(frame.loc[frame["label"] == "T2_reference_dx1_r2", "K2"].iloc[0])
    fine_b0 = float(frame.loc[frame["label"] == "B0_fine_dx075", "K2"].iloc[0])
    legacy_t2 = float(frame.loc[frame["label"] == "T2_fine_legacy_r2", "K2"].iloc[0])
    corrected_t2 = float(frame.loc[frame["label"] == "T2_fine_corrected_r3", "K2"].iloc[0])
    ref_ratio = ref_t2 / ref_b0
    legacy_ratio = legacy_t2 / fine_b0
    corrected_ratio = corrected_t2 / fine_b0

    frame["K2_over_matching_B0"] = np.nan
    frame.loc[frame["label"] == "T2_reference_dx1_r2", "K2_over_matching_B0"] = ref_ratio
    frame.loc[frame["label"] == "T2_fine_legacy_r2", "K2_over_matching_B0"] = legacy_ratio
    frame.loc[frame["label"] == "T2_fine_corrected_r3", "K2_over_matching_B0"] = corrected_ratio

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index=False)

    legacy_error = abs(legacy_ratio / ref_ratio - 1.0)
    corrected_error = abs(corrected_ratio / ref_ratio - 1.0)
    improved = corrected_error < legacy_error
    report = output.with_suffix(".md")
    report.write_text(
        "# Minimal TJ gate-radius test\n\n"
        f"Reference T2/B0 K2 ratio: {ref_ratio:.6f}\n\n"
        f"Fine-grid legacy r=2 T2/B0 ratio: {legacy_ratio:.6f} "
        f"(relative deviation {legacy_error:.2%})\n\n"
        f"Fine-grid corrected r=3 T2/B0 ratio: {corrected_ratio:.6f} "
        f"(relative deviation {corrected_error:.2%})\n\n"
        f"Gate-footprint hypothesis improved convergence: **{improved}**\n\n"
        "The r=3 case is a deliberately minimal pixel-radius test, not the final API. "
        "If it materially improves convergence, replace tj_correlation_radius with a "
        "physical-length parameter in the production model and regression-test the "
        "discretized footprint across dx.\n",
        encoding="utf-8",
    )
    print(frame.to_string(index=False))
    print(report)


if __name__ == "__main__":
    main()
