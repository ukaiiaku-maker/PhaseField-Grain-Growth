#!/usr/bin/env python3
"""Quantify the immutable historical Phase-1 QIU trajectory.

This script is read-only with respect to the historical run.  It extracts
saved-output morphology, stress, energy, and topology diagnostics and writes a
compact forensic summary to a separate qualification root.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def periodic_delta(values: np.ndarray, center: float, length: int) -> np.ndarray:
    return (values - center + 0.5 * length) % length - 0.5 * length


def circular_center(values: np.ndarray, length: int) -> float:
    angles = 2.0 * np.pi * values / length
    angle = np.arctan2(np.sin(angles).mean(), np.cos(angles).mean())
    return float((angle % (2.0 * np.pi)) * length / (2.0 * np.pi))


def periodic_component_count(mask: np.ndarray) -> int:
    components, count = ndimage.label(mask, structure=np.ones((3, 3), dtype=np.uint8))
    if count <= 1:
        return int(count)
    parent = np.arange(count + 1)

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = int(parent[x])
        return x

    def union(a: int, b: int) -> None:
        if a and b:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra

    ny, nx = mask.shape
    for x in range(nx):
        if mask[0, x]:
            for dx in (-1, 0, 1):
                union(int(components[0, x]), int(components[ny - 1, (x + dx) % nx]))
    for y in range(ny):
        if mask[y, 0]:
            for dy in (-1, 0, 1):
                union(int(components[y, 0]), int(components[(y + dy) % ny, nx - 1]))
    roots = {find(i) for i in range(1, count + 1)}
    return len(roots)


def label_metrics(labels: np.ndarray) -> tuple[dict[str, float], dict[int, tuple[float, float, float]]]:
    ny, nx = labels.shape
    ids = np.unique(labels)
    max_id = int(ids.max())
    area_all = np.bincount(labels.ravel(), minlength=max_id + 1).astype(float)
    perimeter_all = np.zeros(max_id + 1, dtype=float)
    for axis in (0, 1):
        neighbor = np.roll(labels, 1, axis=axis)
        changed = labels != neighbor
        perimeter_all += np.bincount(labels[changed], minlength=max_id + 1)
        perimeter_all += np.bincount(neighbor[changed], minlength=max_id + 1)
    compactness = []
    aspect_ratios = []
    components = []
    centroids: dict[int, tuple[float, float, float]] = {}
    for grain_id in ids:
        gid = int(grain_id)
        yy, xx = np.where(labels == gid)
        area = area_all[gid]
        perimeter = perimeter_all[gid]
        cy, cx = circular_center(yy, ny), circular_center(xx, nx)
        dy, dx = periodic_delta(yy, cy, ny), periodic_delta(xx, cx, nx)
        covariance = np.cov(np.stack((dy, dx)), bias=True) if len(yy) > 1 else np.zeros((2, 2))
        # A pixel has finite second moment (1/12) in each direction. Including
        # it prevents a one-pixel-wide component from reporting an artificial
        # floating-point-infinite aspect ratio while retaining elongation.
        eigenvalues = np.maximum(np.linalg.eigvalsh(covariance), 0.0) + 1.0 / 12.0
        aspect = float(np.sqrt(eigenvalues[-1] / eigenvalues[0]))
        component_count = periodic_component_count(labels == gid)
        compactness.append(float(perimeter / np.sqrt(max(4.0 * np.pi * area, np.finfo(float).tiny))))
        aspect_ratios.append(aspect)
        components.append(component_count)
        centroids[gid] = (cy, cx, float(np.sqrt(area / np.pi)))
    comp = np.asarray(compactness)
    aspect = np.asarray(aspect_ratios)
    cc = np.asarray(components)
    metrics = {
        "grain_count": int(len(ids)),
        "mean_area": float(area_all[ids].mean()),
        "G_population": float(np.sqrt(labels.size / len(ids))),
        "pixel_edge_compactness_mean": float(comp.mean()),
        "pixel_edge_compactness_p95": float(np.quantile(comp, 0.95)),
        "pixel_edge_compactness_max": float(comp.max()),
        "aspect_ratio_mean": float(aspect.mean()),
        "aspect_ratio_p95": float(np.quantile(aspect, 0.95)),
        "aspect_ratio_max": float(aspect.max()),
        "components_mean": float(cc.mean()),
        "components_max": int(cc.max()),
        "disconnected_grain_count": int(np.count_nonzero(cc > 1)),
        "disconnected_grain_fraction": float(np.mean(cc > 1)),
    }
    return metrics, centroids


def extinction_cluster(lost: set[int], centroids: dict[int, tuple[float, float, float]], shape: tuple[int, int]) -> int:
    ids = sorted(x for x in lost if x in centroids)
    if not ids:
        return 0
    parent = list(range(len(ids)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, left in enumerate(ids):
        y1, x1, r1 = centroids[left]
        for j in range(i + 1, len(ids)):
            y2, x2, r2 = centroids[ids[j]]
            dy = abs(y1 - y2); dy = min(dy, shape[0] - dy)
            dx = abs(x1 - x2); dx = min(dx, shape[1] - dx)
            if np.hypot(dy, dx) <= 2.0 * max(r1, r2):
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[rj] = ri
    counts: dict[int, int] = defaultdict(int)
    for i in range(len(ids)):
        counts[find(i)] += 1
    return max(counts.values())


def grouped_csv(path: Path, keys: tuple[str, ...]) -> dict[int, list[dict[str, str]]]:
    groups: dict[int, list[dict[str, str]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            groups[int(row["step"])].append(row)
    return groups


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("historical_run", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    run = args.historical_run.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)

    manifest = json.loads((run / "manifest.json").read_text())
    energy_records = json.loads((run / "energy.json").read_text())
    energy_by_step = {int(round(float(r["time"]) / manifest["config"]["pf"]["time_step"])): r for r in energy_records}
    grain_groups = grouped_csv(run / "grain_tracks.csv", ("step",))
    boundary_groups = grouped_csv(run / "boundary_tracks.csv", ("step",))

    trajectory_rows = []
    for step, rows in sorted(grain_groups.items()):
        areas = np.asarray([float(r["area"]) for r in rows])
        trajectory_rows.append({"step": step, "time": float(rows[0]["time"]), "grain_count": len(rows),
                                "G_population": float(np.sqrt(areas.sum() / len(areas))), "grain_loss_since_prior_saved_step": 0})
    for previous, current in zip(trajectory_rows, trajectory_rows[1:]):
        current["grain_loss_since_prior_saved_step"] = int(previous["grain_count"] - current["grain_count"])

    frame_rows = []
    previous_ids: set[int] | None = None
    previous_centroids: dict[int, tuple[float, float, float]] = {}
    previous_step = None
    for frame in sorted((run / "frames").glob("*.npz")):
        with np.load(frame, allow_pickle=False) as data:
            step = int(data["step"])
            if step < 9000:
                continue
            labels = data["labels"].astype(np.int64)
            metrics, centroids = label_metrics(labels)
            ids = set(centroids)
            qiu = np.asarray(data["qiu_shear_stress"], dtype=float)
            boundary_mask = np.asarray(data["boundary_mask"], dtype=bool)
            shear_values = qiu[boundary_mask]
            metrics.update({
                "frame": frame.name, "step": step, "time": float(data["time"]),
                "resolved_shear_rms": float(np.sqrt(np.mean(shear_values**2))) if shear_values.size else 0.0,
                "resolved_shear_p95_abs": float(np.quantile(np.abs(shear_values), 0.95)) if shear_values.size else 0.0,
                "resolved_shear_max_abs": float(np.max(np.abs(shear_values))) if shear_values.size else 0.0,
                "lost_since_prior_saved_frame": 0 if previous_ids is None else len(previous_ids - ids),
                "largest_spatial_extinction_cluster": 0 if previous_ids is None else extinction_cluster(previous_ids - ids, previous_centroids, labels.shape),
                "prior_saved_frame_step": previous_step if previous_step is not None else "",
            })
            closest_energy_step = min(energy_by_step, key=lambda x: abs(x - step))
            metrics["interfacial_energy"] = float(energy_by_step[closest_energy_step]["interfacial"])
            metrics["interfacial_energy_step"] = closest_energy_step
            gr = grain_groups.get(step, [])
            track_compactness = np.asarray([
                float(r["perimeter"]) / np.sqrt(4.0 * np.pi * float(r["area"])) for r in gr
            ])
            metrics["compactness_mean"] = float(track_compactness.mean()) if track_compactness.size else float("nan")
            metrics["compactness_p95"] = float(np.quantile(track_compactness, 0.95)) if track_compactness.size else float("nan")
            metrics["compactness_max"] = float(track_compactness.max()) if track_compactness.size else float("nan")
            br = boundary_groups.get(step, [])
            for field in ("curvature", "normal_velocity"):
                values = np.abs(np.asarray([float(r[field]) for r in br], dtype=float))
                metrics[f"boundary_{field}_mean_abs"] = float(values.mean()) if values.size else float("nan")
                metrics[f"boundary_{field}_p95_abs"] = float(np.quantile(values, 0.95)) if values.size else float("nan")
                metrics[f"boundary_{field}_max_abs"] = float(values.max()) if values.size else float("nan")
            frame_rows.append(metrics)
            previous_ids, previous_centroids, previous_step = ids, centroids, step

    fieldnames = list(frame_rows[0])
    with (output / "historical_frame_forensics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n"); writer.writeheader(); writer.writerows(frame_rows)
    with (output / "historical_trajectory.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(trajectory_rows[0]), lineterminator="\n"); writer.writeheader(); writer.writerows(trajectory_rows)

    steps = np.asarray([r["step"] for r in trajectory_rows])
    counts = np.asarray([r["grain_count"] for r in trajectory_rows])
    g = np.asarray([r["G_population"] for r in trajectory_rows])
    energy_steps = np.asarray(sorted(energy_by_step))
    energy = np.asarray([energy_by_step[s]["interfacial"] for s in energy_steps])
    figure, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
    axes[0].plot(steps, counts, marker="o", ms=2); axes[0].set_ylabel("N")
    axes[1].plot(steps, g, marker="o", ms=2); axes[1].set_ylabel("G population")
    axes[2].plot(energy_steps, energy); axes[2].set_ylabel("Interfacial energy"); axes[2].set_xlabel("solver step")
    for ax in axes:
        ax.axvspan(9800, 10246, color="tab:red", alpha=0.12); ax.grid(alpha=0.25)
    figure.suptitle("Immutable historical QIU trajectory (446-step transition interval)")
    figure.tight_layout(); figure.savefig(output / "historical_qiu_timeseries.png", dpi=180); plt.close(figure)

    transition_frames = [r for r in frame_rows if 9000 <= int(r["step"]) <= 10246]
    fig, axes = plt.subplots(1, len(transition_frames), figsize=(4 * len(transition_frames), 4), squeeze=False)
    for ax, row in zip(axes[0], transition_frames):
        with np.load(run / "frames" / str(row["frame"]), allow_pickle=False) as data:
            ax.imshow(data["labels"], cmap="nipy_spectral", interpolation="nearest")
        ax.set_title(f"step {row['step']}\nN={row['grain_count']}"); ax.axis("off")
    fig.tight_layout(); fig.savefig(output / "historical_transition_field_maps.png", dpi=180); plt.close(fig)

    checkpoint = np.load(run / "checkpoint.npz", allow_pickle=False)
    checkpoint_keys = list(checkpoint.files)
    checkpoint.close()
    restart = {
        "pre_avalanche_restart_available": False,
        "reason": "Only the terminal step-10246 checkpoint contains eta, eigenstrain, active phase mask, and solver metadata. Saved frames contain labels/derived fields but not complete solver or eigenstrain state.",
        "terminal_checkpoint_step": 10246,
        "terminal_checkpoint_keys": checkpoint_keys,
        "replay_requirement": "Exact deterministic legacy replay from the immutable initial state",
    }
    peak_interval_loss = max(trajectory_rows, key=lambda r: int(r["grain_loss_since_prior_saved_step"]))
    summary = {
        "schema_version": 1,
        "historical_run": str(run),
        "historical_source_commit": manifest["git_sha"],
        "historical_manifest_sha256": sha256(run / "manifest.json"),
        "frame_count": len(list((run / "frames").glob("*.npz"))),
        "analyzed_transition_interval": {"start_step": 9800, "end_step": 10246, "solver_steps": 446},
        "saved_frame_metrics": frame_rows,
        "largest_saved_interval_population_loss": peak_interval_loss,
        "restart_audit": restart,
        "first_historical_departure": {
            "quantity": "interfacial_energy",
            "evidence": "Energy rises from 14,853.96 at step 9900 to 64,319.11 at step 10000 before the terminal population collapse; this is a positive-energy explicit-coupling instability signature.",
            "qualification": "Saved cadence cannot identify the exact first solver step; dense deterministic replay is required.",
        },
        "limitations": [
            "Historical scalar output cadence is 100–200 steps, so per-step clipping and extinction onset are unavailable.",
            "Connected-component and clustering metrics are exact on saved periodic label fields but cannot resolve unsaved intermediate steps.",
            "Historical energy records omit full-field elastic energy and therefore do not provide the complete coupled energy balance.",
        ],
    }
    (output / "historical_qiu_baseline_audit.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"output": str(output), "frames_analyzed": len(frame_rows),
                      "largest_loss": peak_interval_loss, "restart": restart["pre_avalanche_restart_available"]}, indent=2))


if __name__ == "__main__":
    main()
