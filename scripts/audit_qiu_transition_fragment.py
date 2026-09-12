#!/usr/bin/env python3
"""Audit and causally reconstruct a bounded legacy-Qiu transition fragment."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyarrow.dataset as ds

from analyze_qiu_full_field_qualification import (
    pre_extinction_precursor,
    transition_timing,
)
from grain_growth_pf.entities.arclength_tracker import ArclengthEntityTracker
from grain_growth_pf.io.checkpoints import atomic_write_text


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def field_path(run: Path, step: int) -> Path:
    paths = sorted((run / "diagnostic_fields").glob(f"step-{step:07d}-*.npz"))
    if not paths:
        raise FileNotFoundError(f"no diagnostic field for step {step}")
    return paths[0]


def load_field(run: Path, step: int) -> dict[str, np.ndarray]:
    path = field_path(run, step)
    with np.load(path) as data:
        return {
            key: np.asarray(data[key]).copy()
            for key in ("eta", "labels", "source_increment")
        }


def domain_measurements(
    phase: np.ndarray,
    previous_phase: np.ndarray,
    partner_phase: np.ndarray,
    points: np.ndarray,
    elapsed: float,
    dx: float,
    *,
    periodic: bool,
) -> dict[str, object]:
    """Expose the two independent gradient stencils used by legacy velocity."""
    coordinates = np.asarray(points, dtype=float).astype(int)
    height, width = phase.shape
    y, x = coordinates[:, 0], coordinates[:, 1]

    def index(values: np.ndarray, size: int) -> np.ndarray:
        return values % size if periodic else np.clip(values, 0, size - 1)

    yc, xc = index(y, height), index(x, width)
    radius = 2
    yp, ym = index(y + radius, height), index(y - radius, height)
    xp, xm = index(x + radius, width), index(x - radius, width)
    span = 2.0 * radius * dx
    gradient_radius2 = np.hypot(
        (phase[yp, xc] - phase[ym, xc]) / span,
        (phase[yc, xp] - phase[yc, xm]) / span,
    )
    yp1, ym1 = index(y + 1, height), index(y - 1, height)
    xp1, xm1 = index(x + 1, width), index(x - 1, width)
    gradient_radius1 = np.hypot(
        (phase[yp1, xc] - phase[ym1, xc]) / (2.0 * dx),
        (phase[yc, xp1] - phase[yc, xm1]) / (2.0 * dx),
    )
    delta_eta = phase[yc, xc] - previous_phase[yc, xc]
    valid = gradient_radius2 > 1e-10
    valid &= (phase[yc, xc] > 0.05) & (phase[yc, xc] < 0.95)
    valid &= phase[yc, xc] + partner_phase[yc, xc] > 0.90
    point_velocity = delta_eta / elapsed / np.maximum(gradient_radius1, 1e-14)
    rows = []
    for index_value in range(len(coordinates)):
        rows.append({
            "point": [int(y[index_value]), int(x[index_value])],
            "valid": bool(valid[index_value]),
            "eta": float(phase[yc[index_value], xc[index_value]]),
            "partner_eta": float(partner_phase[yc[index_value], xc[index_value]]),
            "delta_eta": float(delta_eta[index_value]),
            "gradient_radius2": float(gradient_radius2[index_value]),
            "gradient_radius1": float(gradient_radius1[index_value]),
            "point_velocity": float(point_velocity[index_value]),
        })
    return {
        "points": rows,
        "valid_point_count": int(np.count_nonzero(valid)),
        "mean_velocity": float(np.mean(point_velocity[valid])) if np.any(valid) else 0.0,
        "minimum_radius1_gradient_at_valid_point": (
            float(np.min(gradient_radius1[valid])) if np.any(valid) else None
        ),
        "minimum_radius2_gradient_at_valid_point": (
            float(np.min(gradient_radius2[valid])) if np.any(valid) else None
        ),
    }


def reconstruct_source(
    current: dict[str, np.ndarray],
    previous: dict[str, np.ndarray],
    *,
    elapsed: float,
    dx: float,
    domain_length: float,
    beta: float,
) -> tuple[np.ndarray, dict[str, object]]:
    # Import here so the public stencil audit above remains independently testable.
    from grain_growth_pf.pf.kinematics import interface_kinematics

    eta = current["eta"]
    labels = current["labels"]
    tracker = ArclengthEntityTracker(
        np.zeros(eta.shape[0]), dx=dx, domain_length=domain_length, periodic=True
    )
    snapshot = tracker.update(labels)
    predicted = np.zeros_like(current["source_increment"])
    largest: tuple[float, str, object, float, float, tuple[float, float], np.ndarray] | None = None
    for key, segment in snapshot.boundaries.items():
        _, velocity, normal_tuple = interface_kinematics(
            eta[segment.grain_i], previous["eta"][segment.grain_i],
            segment.points, elapsed, dx, periodic=True,
            partner_phase=eta[segment.grain_j],
        )
        displacement = float(velocity * elapsed)
        normal = np.asarray(normal_tuple, dtype=float)
        tangent = np.asarray((-normal[1], normal[0]))
        vector = beta * displacement * tangent
        strain = 0.5 * (np.outer(vector, normal) + np.outer(normal, vector))
        position = tuple(segment.points[len(segment.points) // 2].astype(int))
        predicted[:, :, position[0] % labels.shape[0], position[1] % labels.shape[1]] += strain
        norm = float(np.linalg.norm(strain))
        if largest is None or norm > largest[0]:
            largest = (norm, key, segment, float(velocity), displacement, normal_tuple, strain)
    assert largest is not None
    norm, key, segment, velocity, displacement, normal, strain = largest
    measurements = domain_measurements(
        eta[segment.grain_i], previous["eta"][segment.grain_i],
        eta[segment.grain_j], segment.points, elapsed, dx, periodic=True,
    )
    position = tuple(segment.points[len(segment.points) // 2].astype(int))
    return predicted, {
        "boundary_entity_id": key,
        "grain_i": int(segment.grain_i),
        "grain_j": int(segment.grain_j),
        "domain_length_pixels": int(len(segment.points)),
        "deposit_position": [int(position[0]), int(position[1])],
        "velocity": velocity,
        "normal_displacement": displacement,
        "normal": [float(normal[0]), float(normal[1])],
        "source_tensor": strain.tolist(),
        "source_tensor_l2": norm,
        "stencil_measurements": measurements,
    }


def relative_l2(residual: np.ndarray, reference: np.ndarray) -> float:
    return float(np.linalg.norm(residual) / max(np.linalg.norm(reference), np.finfo(float).tiny))


def causal_interpretation(*, accumulator_reset: bool) -> str:
    source_state = (
        "The closure reset the diagnostic source accumulator, so the current "
        "field can be compared directly with the reconstructed step source."
        if accumulator_reset else
        "The closure bypassed begin_source_step, so the diagnostic source "
        "accumulator retains prior increments."
    )
    return (
        "A radius-two-valid interface sample has a zero radius-one velocity "
        "gradient. The 1e-14 denominator floor converts finite phase change into "
        "an enormous local velocity and eigenstrain source before stress feedback "
        f"or morphology failure. {source_state}"
    )


def audit(args: argparse.Namespace) -> dict[str, object]:
    run = args.run.resolve()
    scalar_paths = sorted((run / "per_step_diagnostics.parquet").glob("part-*.parquet"))
    scalar = ds.dataset(run / "per_step_diagnostics.parquet", format="parquet").to_table().to_pandas()
    scalar = scalar.sort_values(["step", "time"])
    steps = np.asarray(scalar["step"], dtype=int)
    unique_steps = np.unique(steps)
    expected_steps = np.arange(args.first_step, args.last_step + 1)
    current = load_field(run, args.transition_step)
    previous = load_field(run, args.transition_step - 1)
    predicted, maximum_event = reconstruct_source(
        current, previous, elapsed=args.elapsed, dx=args.dx,
        domain_length=args.domain_length, beta=args.beta,
    )
    observed_delta = current["source_increment"] - previous["source_increment"]
    delta_residual = observed_delta - predicted
    current_residual = current["source_increment"] - predicted
    prior_source_norm = float(np.linalg.norm(previous["source_increment"]))
    retained_residual_norm = float(np.linalg.norm(current_residual))
    accumulator_reset = retained_residual_norm <= 1e-8
    field_paths = sorted((run / "diagnostic_fields").glob("step-*.npz"))
    field_steps = sorted({int(path.name.split("-")[1]) for path in field_paths})
    checkpoint_json = run / "checkpoint.json"
    checkpoint_npz = run / "checkpoint.npz"
    checkpoint = json.loads(checkpoint_json.read_text())
    with np.load(checkpoint_npz) as data:
        checkpoint_arrays = sorted(data.files)

    selected_rows = {}
    timing = transition_timing(scalar)
    for role, step in (
        ("pre_transition", args.transition_step - 1),
        ("transition", args.transition_step),
        ("feedback", args.transition_step + 1),
        ("first_disconnection", timing["first_disconnected_grain_step"]),
        ("fragment_terminal", args.last_step),
    ):
        if step is None:
            continue
        row = scalar.loc[scalar["step"] == step]
        if row.empty:
            continue
        record = row.iloc[-1]
        selected_rows[role] = {
            key: (int(record[key]) if key in {
                "step", "grain_count", "newly_extinct_phases",
                "largest_100_step_population_loss", "disconnected_grain_count",
            } else float(record[key]))
            for key in (
                "step", "time", "grain_count", "newly_extinct_phases",
                "largest_100_step_population_loss", "source_increment_l2",
                "stress_linf", "eigenstrain_linf", "interfacial_energy",
                "elastic_energy", "total_energy", "clipped_fraction",
                "compactness_mean", "compactness_max", "aspect_ratio_max",
                "disconnected_grain_count", "max_external_pair_driving_difference",
            )
        }

    archive = args.archive.resolve() if args.archive else None
    tracker_manifest = args.tracker_manifest.resolve() if args.tracker_manifest else None
    result: dict[str, object] = {
        "schema_version": 1,
        "scientific_status": "bounded_nonterminal_transition_fragment_audited",
        "run": str(run),
        "archive": ({
            "path": str(archive), "sha256": sha256(archive),
            "size_bytes": archive.stat().st_size,
        } if archive else None),
        "checkpoint": {
            "restart_capable_through_step": int(checkpoint["step_number"]),
            "json_path": str(checkpoint_json), "json_sha256": sha256(checkpoint_json),
            "npz_path": str(checkpoint_npz), "npz_sha256": sha256(checkpoint_npz),
            "npz_arrays": checkpoint_arrays,
        },
        "scalar_stream": {
            "part_count": len(scalar_paths), "row_count": int(len(scalar)),
            "unique_step_count": int(len(unique_steps)),
            "step_first": int(unique_steps[0]), "step_last": int(unique_steps[-1]),
            "duplicate_step_count": int(len(steps) - len(unique_steps)),
            "missing_expected_steps": sorted(set(expected_steps) - set(unique_steps)),
        },
        "field_stream": {
            "file_count": len(field_paths), "unique_step_count": len(field_steps),
            "step_first": field_steps[0], "step_last": field_steps[-1],
            "dense_window_first": args.dense_field_start,
            "dense_window_missing_steps": sorted(
                set(range(args.dense_field_start, args.last_step + 1)) - set(field_steps)
            ),
        },
        "tracker_window_manifest": ({
            "path": str(tracker_manifest), "sha256": sha256(tracker_manifest),
            "content": json.loads(tracker_manifest.read_text()),
        } if tracker_manifest else None),
        "transition_timing": timing,
        "pre_extinction_precursor": pre_extinction_precursor(scalar),
        "selected_rows": selected_rows,
        "causal_reconstruction": {
            "transition_step": args.transition_step,
            "previous_field_sha256": sha256(field_path(run, args.transition_step - 1)),
            "transition_field_sha256": sha256(field_path(run, args.transition_step)),
            "observed_source_delta_l2": float(np.linalg.norm(observed_delta)),
            "predicted_source_l2": float(np.linalg.norm(predicted)),
            "previous_accumulator_l2": prior_source_norm,
            "current_minus_predicted_l2": retained_residual_norm,
            "observed_delta_minus_predicted_relative_l2": relative_l2(delta_residual, predicted),
            "current_accumulator_minus_predicted_relative_l2": relative_l2(current_residual, predicted),
            "source_accumulator_was_reset_by_closure": bool(accumulator_reset),
            "largest_reconstructed_event": maximum_event,
            "causal_interpretation": causal_interpretation(
                accumulator_reset=accumulator_reset
            ),
        },
    }
    return result


def plot_history(run: Path, output: Path, transition_step: int) -> None:
    frame = ds.dataset(run / "per_step_diagnostics.parquet", format="parquet").to_table().to_pandas()
    frame = frame.sort_values("step").drop_duplicates("step", keep="last")
    figure, axes = plt.subplots(3, 1, figsize=(9.0, 9.0), sharex=True)
    for field in ("source_increment_l2", "stress_linf", "eigenstrain_linf"):
        axes[0].semilogy(frame["step"], np.maximum(np.abs(frame[field]), 1e-30), label=field)
    axes[0].set_ylabel("source / field norm"); axes[0].legend(fontsize=8)
    for field in ("interfacial_energy", "elastic_energy", "total_energy"):
        axes[1].semilogy(frame["step"], np.maximum(np.abs(frame[field]), 1e-30), label=field)
    axes[1].set_ylabel("energy magnitude"); axes[1].legend(fontsize=8)
    axes[2].plot(frame["step"], frame["grain_count"], label="grain_count")
    axes[2].plot(frame["step"], frame["disconnected_grain_count"], label="disconnected_grains")
    axes[2].plot(frame["step"], frame["compactness_mean"], label="mean_compactness")
    axes[2].set_ylabel("morphology"); axes[2].set_xlabel("step"); axes[2].legend(fontsize=8)
    for axis in axes:
        axis.axvline(transition_step, color="black", linestyle="--", linewidth=1.0)
        axis.grid(alpha=0.2)
    figure.suptitle("Legacy Qiu replay: source singularity precedes morphology cascade")
    figure.tight_layout(); figure.savefig(output, dpi=180); plt.close(figure)


def plot_local_cause(run: Path, output: Path, transition_step: int) -> None:
    previous = load_field(run, transition_step - 1)
    current = load_field(run, transition_step)
    source_delta = current["source_increment"] - previous["source_increment"]
    source_magnitude = np.sqrt(np.sum(source_delta * source_delta, axis=(0, 1)))
    center = np.unravel_index(int(np.argmax(source_magnitude)), source_magnitude.shape)
    center_pair = (int(center[0]), int(center[1]))
    radius = 14
    yy = (np.arange(center[0] - radius, center[0] + radius + 1) % source_magnitude.shape[0])
    xx = (np.arange(center[1] - radius, center[1] + radius + 1) % source_magnitude.shape[1])
    labels = current["labels"][np.ix_(yy, xx)]
    phase = current["eta"][343][np.ix_(yy, xx)]
    delta_eta = (current["eta"][343] - previous["eta"][343])[np.ix_(yy, xx)]
    source = source_magnitude[np.ix_(yy, xx)]
    figure, axes = plt.subplots(1, 4, figsize=(13.0, 3.5), constrained_layout=True)
    panels = (
        (labels, "grain labels", "tab20"),
        (phase, "$\\eta_{343}$", "viridis"),
        (delta_eta, "$\\Delta\\eta_{343}$", "coolwarm"),
        (np.log10(np.maximum(source, 1e-30)), "$\\log_{10}||\\Delta\\epsilon^*||$", "magma"),
    )
    for axis, (values, title, cmap) in zip(axes, panels):
        image = axis.imshow(values, origin="lower", cmap=cmap)
        axis.set_title(title); axis.set_xticks([]); axis.set_yticks([])
        figure.colorbar(image, ax=axis, fraction=0.046)
    axes[-1].scatter([radius], [radius], facecolors="none", edgecolors="cyan", s=80, linewidths=1.5)
    figure.suptitle(
        f"Step {transition_step} localized source singularity near pixel {center_pair}"
    )
    figure.savefig(output, dpi=180); plt.close(figure)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("run", type=Path)
    result.add_argument("output", type=Path)
    result.add_argument("--archive", type=Path)
    result.add_argument("--tracker-manifest", type=Path)
    result.add_argument("--first-step", type=int, required=True)
    result.add_argument("--last-step", type=int, required=True)
    result.add_argument("--dense-field-start", type=int, required=True)
    result.add_argument("--transition-step", type=int, required=True)
    result.add_argument("--elapsed", type=float, default=0.04)
    result.add_argument("--dx", type=float, default=1.0)
    result.add_argument("--domain-length", type=float, default=12.0)
    result.add_argument("--beta", type=float, default=0.35)
    return result


def main() -> None:
    args = parser().parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    result = audit(args)
    audit_path = args.output / "transition_fragment_audit.json"
    atomic_write_text(audit_path, json.dumps(result, indent=2) + "\n")
    plot_path = args.output / "transition_precursor_history.png"
    plot_history(args.run.resolve(), plot_path, args.transition_step)
    plot_local_cause(
        args.run.resolve(), args.output / "transition_local_cause.png",
        args.transition_step,
    )
    checksums = {
        path.name: sha256(path) for path in sorted(args.output.iterdir()) if path.is_file()
    }
    atomic_write_text(args.output / "checksums.json", json.dumps(checksums, indent=2) + "\n")
    print(audit_path.resolve())


if __name__ == "__main__":
    main()
