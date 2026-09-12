#!/usr/bin/env python3
"""Audit per-boundary/source-work closure across a Qiu diagnostic stream."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

from analyze_qiu_full_field_qualification import json_safe, sha256, transition_timing
from grain_growth_pf.io.checkpoints import atomic_write_text


def maximum_relative_difference(first: pd.Series, second: pd.Series) -> float:
    left = np.asarray(first, dtype=float)
    right = np.asarray(second, dtype=float)
    scale = np.maximum(np.maximum(np.abs(left), np.abs(right)), np.finfo(float).tiny)
    return float(np.max(np.abs(left - right) / scale))


def reconcile_streams(
    scalar: pd.DataFrame, boundary: pd.DataFrame, *, transition_step: int
) -> tuple[pd.DataFrame, dict[str, object]]:
    aggregates = boundary.groupby("step", as_index=False).agg(
        boundary_predicted_work=("predicted_elastic_work", "sum"),
        boundary_integrated_sweep=("integrated_sweep", "sum"),
    )
    merged = scalar.merge(aggregates, on="step", how="left", validate="one_to_one")
    if merged[["boundary_predicted_work", "boundary_integrated_sweep"]].isna().any().any():
        raise RuntimeError("one or more scalar steps lack per-boundary diagnostics")

    work_difference = np.abs(
        merged["predicted_source_work"] - merged["boundary_predicted_work"]
    )
    sweep_difference = np.abs(
        merged["integrated_boundary_sweep"] - merged["boundary_integrated_sweep"]
    )
    pre = merged.loc[merged["step"] < transition_step]
    return merged, {
        "all_steps": {
            "maximum_abs_predicted_work_difference": float(work_difference.max()),
            "maximum_relative_predicted_work_difference": maximum_relative_difference(
                merged["predicted_source_work"], merged["boundary_predicted_work"]
            ),
            "maximum_abs_integrated_sweep_difference": float(sweep_difference.max()),
            "maximum_relative_integrated_sweep_difference": maximum_relative_difference(
                merged["integrated_boundary_sweep"], merged["boundary_integrated_sweep"]
            ),
        },
        "before_transition": {
            "step_count": int(len(pre)),
            "maximum_abs_predicted_work_difference": float(
                np.abs(pre["predicted_source_work"] - pre["boundary_predicted_work"]).max()
            ),
            "maximum_relative_predicted_work_difference": maximum_relative_difference(
                pre["predicted_source_work"], pre["boundary_predicted_work"]
            ),
            "maximum_abs_integrated_sweep_difference": float(
                np.abs(pre["integrated_boundary_sweep"] - pre["boundary_integrated_sweep"]).max()
            ),
            "maximum_relative_integrated_sweep_difference": maximum_relative_difference(
                pre["integrated_boundary_sweep"], pre["boundary_integrated_sweep"]
            ),
        },
    }


def _row(frame: pd.DataFrame, step: int) -> dict[str, object]:
    selected = frame.loc[frame["step"] == step]
    if selected.empty:
        raise RuntimeError(f"missing selected step {step}")
    fields = (
        "step", "time", "grain_count", "newly_extinct_phases",
        "source_increment_l2", "source_increment_linf",
        "source_elastic_energy_change", "predicted_source_work",
        "boundary_predicted_work", "source_work_error", "stress_linf",
        "eigenstrain_linf", "total_energy", "compactness_mean",
        "compactness_max", "aspect_ratio_max", "disconnected_grain_count",
        "largest_100_step_population_loss",
    )
    return {field: selected.iloc[-1][field] for field in fields}


def audit(args: argparse.Namespace) -> dict[str, object]:
    run = args.run.resolve()
    scalar = (
        ds.dataset(run / "per_step_diagnostics.parquet", format="parquet")
        .to_table().to_pandas().sort_values("step").drop_duplicates("step", keep="last")
    )
    boundary = (
        ds.dataset(run / "per_boundary_diagnostics.parquet", format="parquet")
        .to_table().to_pandas()
    )
    steps = np.asarray(scalar["step"], dtype=int)
    expected = np.arange(steps[0], steps[-1] + 1)
    if not np.array_equal(steps, expected):
        raise RuntimeError("scalar step stream is not contiguous")
    merged, reconciliation = reconcile_streams(
        scalar, boundary, transition_step=args.transition_step
    )
    timing = transition_timing(scalar)
    selected_steps = {
        "pre_source_singularity": args.transition_step - 1,
        "source_singularity": args.transition_step,
        "stress_feedback": args.transition_step + 1,
        "first_disconnection": timing["first_disconnected_grain_step"],
        "mean_compactness_gate": timing["first_mean_compactness_above_2p5_step"],
        "max_compactness_gate": timing["first_max_compactness_above_6_step"],
        "population_burst_gate": timing["first_100_step_population_loss_above_10pct_step"],
        "fragment_terminal": int(steps[-1]),
    }
    selected_rows = {
        role: _row(merged, int(step))
        for role, step in selected_steps.items() if step is not None
    }
    transition_boundaries = boundary.loc[boundary["step"] == args.transition_step]
    if transition_boundaries.empty:
        raise RuntimeError("transition step has no per-boundary rows")
    transition_event = transition_boundaries.loc[
        transition_boundaries["source_tensor_norm"].idxmax()
    ].to_dict()
    global_event = boundary.loc[boundary["source_tensor_norm"].idxmax()].to_dict()

    pre = scalar.loc[scalar["step"] < args.transition_step]
    relative_work = pre["source_work_error"].abs() / (
        pre["source_elastic_energy_change"].abs()
        + pre["predicted_source_work"].abs()
        + np.finfo(float).eps
    )
    field_audit = None
    if args.causal_field_audit:
        path = args.causal_field_audit.resolve()
        content = json.loads(path.read_text())
        reconstruction = content["causal_reconstruction"]
        field_event = reconstruction["largest_reconstructed_event"]
        field_audit = {
            "path": str(path),
            "sha256": sha256(path),
            "transition_step": reconstruction["transition_step"],
            "observed_delta_minus_predicted_relative_l2": reconstruction[
                "observed_delta_minus_predicted_relative_l2"
            ],
            "boundary_entity_matches": (
                field_event["boundary_entity_id"] == transition_event["entity_id"]
            ),
            "source_tensor_l2_matches": bool(np.isclose(
                field_event["source_tensor_l2"], transition_event["source_tensor_norm"],
                rtol=1e-14, atol=0.0,
            )),
            "normal_displacement_matches": bool(np.isclose(
                field_event["normal_displacement"], transition_event["normal_displacement"],
                rtol=1e-14, atol=0.0,
            )),
        }

    return {
        "schema_version": 1,
        "scientific_status": "bounded_nonterminal_sourcework_transition_audited",
        "run": str(run),
        "transition_step": args.transition_step,
        "scalar_stream": {
            "row_count": int(len(scalar)), "step_first": int(steps[0]),
            "step_last": int(steps[-1]), "contiguous": True,
        },
        "boundary_stream": {
            "row_count": int(len(boundary)),
            "unique_step_count": int(boundary["step"].nunique()),
            "step_first": int(boundary["step"].min()),
            "step_last": int(boundary["step"].max()),
        },
        "boundary_scalar_reconciliation": reconciliation,
        "transition_timing": timing,
        "selected_rows": selected_rows,
        "largest_transition_boundary": transition_event,
        "largest_fragment_boundary": global_event,
        "pre_transition_source_work_error": {
            "median_absolute": float(pre["source_work_error"].abs().median()),
            "p95_absolute": float(pre["source_work_error"].abs().quantile(0.95)),
            "maximum_absolute": float(pre["source_work_error"].abs().max()),
            "median_relative": float(relative_work.median()),
            "p95_relative": float(relative_work.quantile(0.95)),
            "maximum_relative": float(relative_work.max()),
        },
        "independent_causal_field_audit": field_audit,
        "causal_order_verified": bool(
            args.transition_step < timing["first_disconnected_grain_step"]
            < timing["first_mean_compactness_above_2p5_step"]
            < timing["first_100_step_population_loss_above_10pct_step"]
        ),
    }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("run", type=Path)
    result.add_argument("output", type=Path)
    result.add_argument("--transition-step", type=int, required=True)
    result.add_argument("--causal-field-audit", type=Path)
    return result


def main() -> None:
    args = parser().parse_args()
    result = audit(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_text(
        args.output,
        json.dumps(json_safe(result), indent=2, allow_nan=False) + "\n",
    )
    print(args.output.resolve())


if __name__ == "__main__":
    main()
