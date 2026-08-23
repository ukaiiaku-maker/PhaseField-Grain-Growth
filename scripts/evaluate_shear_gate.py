#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def _ordered_values(frame: pd.DataFrame, metric: str) -> np.ndarray:
    return frame.sort_values("shear_stiffness")[metric].to_numpy(float)


def _middle_deviation(values: np.ndarray) -> float:
    if len(values) != 3 or np.any(~np.isfinite(values)):
        return np.inf
    scale = max(float(np.ptp(values)), 0.1 * float(np.max(np.abs(values))), 1e-14)
    return float(abs(values[1] - 0.5 * (values[0] + values[2])) / scale)


def _same_nonzero_sign(values: list[float], tolerance: float = 1e-12) -> bool:
    signs = [int(np.sign(value)) for value in values if abs(value) > tolerance]
    return len(signs) == len(values) and len(set(signs)) == 1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("analysis")
    parser.add_argument("--roughness", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--conservation-tolerance", type=float, default=1e-9)
    args = parser.parse_args()
    analysis = Path(args.analysis)
    destination = Path(args.output)
    destination.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(analysis / "shear_run_summary.csv")
    windows = pd.read_csv(analysis / "topology_window_metrics.csv")
    causal = pd.read_csv(analysis / "solver_step_causal_nulls.csv")
    roughness = pd.read_csv(args.roughness)
    mandatory = windows[windows["topology_window"] == "N190_to_160"].copy()

    rates = mandatory.set_index("regime")["Rdot"].to_dict()
    gb_differences = [
        float(rates[f"GTSC_GB_Ks{suffix}"] - rates["GTC_GB"])
        for suffix in ("010", "020", "040")
    ]
    gbtj_differences = [
        float(rates[f"GTSC_GBTJ_Ks{suffix}"] - rates[f"GSC_GBTJ_Ks{suffix}"])
        for suffix in ("010", "020", "040")
    ]

    smoothness: dict[str, dict[str, float]] = {}
    monotonic_stress_energy = True
    for family, frame in summary[summary["regime"].str.contains("Ks")].groupby(
        summary["regime"].str.replace(r"_Ks\d+$", "", regex=True)
    ):
        smoothness[str(family)] = {}
        for metric in (
            "RMS_shear_state", "p95_abs_tau_int", "stored_shear_energy_max",
            "p95_abs_tauVtau_over_kBT",
        ):
            values = _ordered_values(frame, metric)
            smoothness[str(family)][metric] = _middle_deviation(values)
        for metric in ("p95_abs_tau_int", "stored_shear_energy_max"):
            values = _ordered_values(frame, metric)
            monotonic_stress_energy &= bool(
                np.all(np.isfinite(values))
                and (
                    np.all(np.diff(values) >= -1e-12)
                    or np.all(np.diff(values) <= 1e-12)
                )
            )

    roughness_checks: dict[str, dict[str, float]] = {}
    bounded_roughness = True
    roughness["family"] = roughness["regime"].str.replace(r"_Ks\d+$", "", regex=True)
    for family, frame in roughness[roughness["regime"].str.contains("Ks")].groupby("family"):
        means = frame.groupby("regime")["isoperimetric_q_mean"].mean()
        cross_scale_ratio = float(means.max() / means.min()) if means.min() > 0 else np.inf
        progressive_ratios = []
        for _, history in frame.groupby("regime"):
            history = history.sort_values("step")
            start = float(history["isoperimetric_q_mean"].iloc[0])
            progressive_ratios.append(
                float(history["isoperimetric_q_mean"].max() / start) if start > 0 else np.inf
            )
        progressive_max = max(progressive_ratios, default=np.inf)
        roughness_checks[str(family)] = {
            "cross_Ks_mean_q_ratio": cross_scale_ratio,
            "max_progressive_q_ratio": progressive_max,
        }
        bounded_roughness &= cross_scale_ratio <= 1.25 and progressive_max <= 1.25

    early_rates = mandatory["Rdot"].to_numpy(float)
    finite_smoothness = [
        value for family in smoothness.values() for value in family.values()
    ]
    checks = {
        "ten_terminal_runs": bool(len(summary) == 10 and (summary["status"] == "completed").all()),
        "mandatory_N190_to_160_available": bool(
            len(mandatory) == 10 and mandatory["available"].eq(True).all()
        ),
        "all_early_growth_rates_positive": bool(
            len(early_rates) == 10 and np.all(np.isfinite(early_rates)) and np.all(early_rates > 0.0)
        ),
        "qualitative_GB_reduced_full_ordering_stable": _same_nonzero_sign(gb_differences),
        "qualitative_GBTJ_reduced_full_ordering_stable": _same_nonzero_sign(gbtj_differences),
        "reference_not_abrupt_threshold": bool(
            finite_smoothness and max(finite_smoothness) <= 0.80
        ),
        "stress_and_energy_scale_continuously": bool(monotonic_stress_energy),
        "no_single_value_pathological_stagnation": bool(
            np.min(early_rates) / np.max(early_rates) >= 0.02
        ),
        "event_traces_present_and_overlap_merged": bool(
            (summary["event_trace_rows"] > 0).all()
            and (summary["event_trace_events"] > 0).all()
            and (summary["duplicate_entity_step_trace_rows"] == 0).all()
        ),
        "event_causal_nulls_complete": bool(
            len(causal) == 10 * 6
            and causal[[
                "circular_shift_null_mean", "block_shuffle_null_mean",
                "same_entity_null_mean", "same_topology_window_null_mean",
            ]].notna().all().all()
        ),
        "conservation_within_tolerance": bool(
            summary["max_abs_conservation_residual"].abs().max()
            <= args.conservation_tolerance
        ),
        "bounded_nonprogressive_morphology": bool(bounded_roughness),
    }
    decision = "GO" if all(checks.values()) else "NO-GO"
    payload: dict[str, Any] = {
        "decision": decision,
        "checks": checks,
        "GB_full_minus_reduced_Rdot": gb_differences,
        "GBTJ_full_minus_reduced_Rdot": gbtj_differences,
        "smoothness_middle_deviation": smoothness,
        "roughness": roughness_checks,
        "conservation_tolerance": args.conservation_tolerance,
        "note": (
            "Automated screen of the stated gate criteria; final scientific gate also "
            "requires visual movie/contact-sheet review."
        ),
    }
    (destination / "shear_gate.json").write_text(json.dumps(payload, indent=2) + "\n")
    lines = [
        f"# 900 K shear-scale decision: {decision}", "",
        "This decision uses the mandatory N=190→160 window and is conditional on visual movie review.", "",
        "## Criteria", "",
    ]
    lines.extend(
        f"- {'PASS' if passed else 'FAIL'} — {name.replace('_', ' ')}"
        for name, passed in checks.items()
    )
    lines.extend(["", "## Numerical diagnostics", "", "```json", json.dumps({
        "GB_full_minus_reduced_Rdot": gb_differences,
        "GBTJ_full_minus_reduced_Rdot": gbtj_differences,
        "smoothness_middle_deviation": smoothness,
        "roughness": roughness_checks,
    }, indent=2), "```", ""])
    (destination / "shear_gate.md").write_text("\n".join(lines))
    print(destination / "shear_gate.md")
    if decision != "GO":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
