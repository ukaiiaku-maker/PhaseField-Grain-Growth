from __future__ import annotations

import json
from collections import deque
from pathlib import Path
from typing import Any

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from numba import njit


@njit(cache=True)
def _component_counts(labels: np.ndarray) -> np.ndarray:
    """Count periodic 8-connected components for every integer label."""
    ny, nx = labels.shape
    size = ny * nx
    parent = np.arange(size)

    def root(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for y in range(ny):
        for x in range(nx):
            here = y * nx + x
            value = labels[y, x]
            for oy, ox in ((0, -1), (-1, -1), (-1, 0), (-1, 1)):
                other_y = (y + oy) % ny
                other_x = (x + ox) % nx
                if labels[other_y, other_x] != value:
                    continue
                left = root(here)
                right = root(other_y * nx + other_x)
                if left != right:
                    parent[right] = left
    maximum = int(np.max(labels))
    counts = np.zeros(maximum + 1, dtype=np.int64)
    for index in range(size):
        if root(index) == index:
            counts[labels[index // nx, index % nx]] += 1
    return counts


def _quantile(values: np.ndarray, q: float) -> float:
    return float(np.quantile(values, q)) if len(values) else float("nan")


def morphology_metrics(labels: np.ndarray, snapshot: Any) -> tuple[dict[str, float | int], dict[int, tuple[float, float, float]]]:
    """Return geometry metrics from the already accepted label field."""
    grains = list(snapshot.grains.values())
    areas = np.asarray([grain.area for grain in grains], dtype=float)
    compactness = np.asarray([
        grain.perimeter / np.sqrt(max(4.0 * np.pi * grain.area, np.finfo(float).tiny))
        for grain in grains
    ])
    components_all = _component_counts(np.asarray(labels, dtype=np.int64))
    components = np.asarray([components_all[grain.grain_id] for grain in grains], dtype=int)

    # Periodic second moments are evaluated about the tracker's circular centroid.
    flat = labels.ravel()
    yy, xx = np.indices(labels.shape)
    aspects: list[float] = []
    for grain in grains:
        mask = flat == grain.grain_id
        y = yy.ravel()[mask]
        x = xx.ravel()[mask]
        dy = (y - grain.centroid[0] + labels.shape[0] / 2) % labels.shape[0] - labels.shape[0] / 2
        dx = (x - grain.centroid[1] + labels.shape[1] / 2) % labels.shape[1] - labels.shape[1] / 2
        covariance = np.cov(np.stack((dy, dx)), bias=True) if len(y) > 1 else np.zeros((2, 2))
        eigenvalues = np.maximum(np.linalg.eigvalsh(covariance), 0.0) + 1.0 / 12.0
        aspects.append(float(np.sqrt(eigenvalues[-1] / eigenvalues[0])))
    aspect = np.asarray(aspects)
    count = len(grains)
    centroids = {
        grain.grain_id: (*grain.centroid, grain.equivalent_radius) for grain in grains
    }
    return {
        "grain_count": count,
        "mean_area": float(np.mean(areas)) if count else float("nan"),
        "area_std": float(np.std(areas)) if count else float("nan"),
        "area_skewness": (
            float(np.mean(((areas - areas.mean()) / areas.std()) ** 3))
            if count and areas.std() > 0 else 0.0
        ),
        "G_population": float(np.sqrt(labels.size / count)) if count else float("nan"),
        "compactness_mean": float(np.mean(compactness)) if count else float("nan"),
        "compactness_p95": _quantile(compactness, 0.95),
        "compactness_max": float(np.max(compactness)) if count else float("nan"),
        "aspect_ratio_mean": float(np.mean(aspect)) if count else float("nan"),
        "aspect_ratio_p95": _quantile(aspect, 0.95),
        "aspect_ratio_max": float(np.max(aspect)) if count else float("nan"),
        "components_mean": float(np.mean(components)) if count else float("nan"),
        "components_max": int(np.max(components)) if count else 0,
        "disconnected_grain_count": int(np.count_nonzero(components > 1)),
    }, centroids


def _extinction_cluster(
    lost: set[int], centroids: dict[int, tuple[float, float, float]], shape: tuple[int, int]
) -> int:
    ids = sorted(lost.intersection(centroids))
    if not ids:
        return 0
    parent = list(range(len(ids)))

    def root(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for left_index, left in enumerate(ids):
        y1, x1, r1 = centroids[left]
        for right_index in range(left_index + 1, len(ids)):
            y2, x2, r2 = centroids[ids[right_index]]
            dy = abs(y1 - y2); dy = min(dy, shape[0] - dy)
            dx = abs(x1 - x2); dx = min(dx, shape[1] - dx)
            if np.hypot(dy, dx) <= 2.0 * max(r1, r2):
                a, b = root(left_index), root(right_index)
                if a != b:
                    parent[b] = a
    counts: dict[int, int] = {}
    for index in range(len(ids)):
        key = root(index)
        counts[key] = counts.get(key, 0) + 1
    return max(counts.values())


class QiuForensicRecorder:
    """Restart-safe, opt-in forensic stream for the legacy full-field path."""

    schema_version = 1

    def __init__(self, output_dir: Path, parameters: dict[str, Any], *, resume: bool = False):
        self.output_dir = Path(output_dir)
        self.scalar_dir = self.output_dir / "per_step_diagnostics.parquet"
        self.boundary_dir = self.output_dir / "per_boundary_diagnostics.parquet"
        self.field_dir = self.output_dir / "diagnostic_fields"
        for path in (self.scalar_dir, self.boundary_dir, self.field_dir):
            path.mkdir(parents=True, exist_ok=True)
        self.flush_rows = max(1, int(parameters.get("qiu_diagnostic_flush_rows", 32)))
        self.field_start = int(parameters.get("qiu_diagnostic_field_start_step", 9000))
        self.field_cadence = max(1, int(parameters.get("qiu_diagnostic_field_cadence", 10)))
        self.clip_spike_fraction = float(parameters.get("qiu_guard_clip_fraction", 0.02))
        self.extinction_spike = int(parameters.get("qiu_guard_extinction_count", 10))
        self.terminate_on_guard = bool(parameters.get("qiu_guard_terminate", True))
        self.scalar_rows: list[dict[str, Any]] = []
        self.boundary_rows: list[dict[str, Any]] = []
        self.scalar_part = self._next_part(self.scalar_dir) if resume else 0
        self.boundary_part = self._next_part(self.boundary_dir) if resume else 0
        self.population_history: deque[tuple[int, int]] = deque()
        self.previous_ids: set[int] = set()
        self.previous_centroids: dict[int, tuple[float, float, float]] = {}
        self.previous_interfacial = float("nan")
        self.before_source_energy = float("nan")
        self.boundary_sweep = 0.0
        self.source_measure = 0.0
        self.predicted_source_work = 0.0
        self.guard_reason: str | None = None

    def initialize(self, simulation: Any, interfacial_energy: float) -> None:
        """Seed differenced diagnostics from the accepted initial/restart state."""
        _, centroids = morphology_metrics(simulation.solver.labels, simulation.snapshot)
        self.previous_ids = set(simulation.snapshot.grains)
        self.previous_centroids = centroids
        self.previous_interfacial = float(interfacial_energy)

    @staticmethod
    def _next_part(directory: Path) -> int:
        parts = sorted(directory.glob("part-*.parquet"))
        return 0 if not parts else int(parts[-1].stem.split("-")[-1]) + 1

    def begin_coupling(self, full_field: Any, grid_spacing: float) -> None:
        self.before_source_energy = full_field.elastic_energy(grid_spacing)
        self.boundary_sweep = 0.0
        self.source_measure = 0.0
        self.predicted_source_work = 0.0

    def record_boundary(
        self, *, step: int, time: float, entity_id: str, grain_i: int, grain_j: int,
        length: float, normal_displacement: float, beta: float,
        resolved_shear: float, source_tensor: np.ndarray,
    ) -> None:
        sweep = float(length * normal_displacement)
        source_integral = float(np.linalg.norm(source_tensor))
        work = float(-resolved_shear * beta * sweep)
        self.boundary_sweep += sweep
        self.source_measure += source_integral
        self.predicted_source_work += work
        self.boundary_rows.append({
            "step": step, "time": time, "entity_id": entity_id,
            "grain_i": grain_i, "grain_j": grain_j, "length": length,
            "normal_displacement": normal_displacement, "beta": beta,
            "resolved_shear": resolved_shear, "integrated_sweep": sweep,
            "source_tensor_norm": source_integral, "predicted_elastic_work": work,
        })

    def record_global_sweep(
        self, *, swept_area: float, integrated_source: np.ndarray,
        predicted_source_work: float,
    ) -> None:
        self.boundary_sweep = float(swept_area)
        self.source_measure = float(np.linalg.norm(integrated_source))
        self.predicted_source_work = float(predicted_source_work)

    def record_step(self, simulation: Any, diag: Any) -> str | None:
        full_field = simulation.full_field
        labels = simulation.solver.labels
        metrics, centroids = morphology_metrics(labels, simulation.snapshot)
        current_ids = set(simulation.snapshot.grains)
        lost = self.previous_ids - current_ids
        cluster = _extinction_cluster(lost, self.previous_centroids, labels.shape)
        count = int(metrics["grain_count"])
        self.population_history.append((diag.step, count))
        while self.population_history and self.population_history[0][0] < diag.step - 100:
            self.population_history.popleft()
        loss_100 = max((old_count - count for _, old_count in self.population_history), default=0)

        eigenstrain = full_field.eigenstrain
        stress = full_field.stress
        source = full_field.source_increment
        abs_stress = np.abs(stress).ravel()
        elastic = full_field.elastic_energy(simulation.config.pf.grid_spacing)
        interfacial = float(diag.interfacial_energy)
        source_delta = elastic - self.before_source_energy
        source_work_error = source_delta - self.predicted_source_work
        capillary_values = np.asarray([
            abs(simulation.config.pf.gb_energy * boundary.curvature)
            for boundary in simulation.snapshot.boundaries.values()
        ])
        elastic_values = np.asarray([
            abs(float(simulation.config.parameters.get("easy_beta", 0.35))
                * simulation._boundary_resolved_shear(
                    simulation.domains[boundary.entity_id], boundary
                ))
            for boundary in simulation.snapshot.boundaries.values()
            if boundary.entity_id in simulation.domains
        ])
        if len(capillary_values) and len(elastic_values):
            length = min(len(capillary_values), len(elastic_values))
            ratios = elastic_values[:length] / (capillary_values[:length] + np.finfo(float).eps)
        else:
            ratios = np.empty(0)
        row: dict[str, Any] = {
            "schema_version": self.schema_version,
            "step": diag.step, "time": diag.time,
            "requested_dt": diag.requested_dt, "used_dt": diag.dt,
            "max_abs_order_parameter_increment": diag.max_abs_order_parameter_increment,
            "max_raw_order_parameter_increment": diag.max_raw_order_parameter_increment,
            "clipped_low_count": diag.clipped_low_count,
            "clipped_high_count": diag.clipped_high_count,
            "trial_value_count": diag.trial_value_count,
            "clipped_fraction": diag.clipped_fraction,
            "total_renormalization_correction": diag.total_renormalization_correction,
            "simplex_residual": diag.max_constraint_error,
            "newly_extinct_phases": diag.newly_extinct_phases,
            "minimum_pre_extinction_maximum": diag.minimum_pre_extinction_maximum,
            "max_external_pair_driving_difference": diag.max_external_pair_driving_difference,
            "eigenstrain_l2": float(np.linalg.norm(eigenstrain)),
            "eigenstrain_linf": float(np.max(np.abs(eigenstrain))),
            "stress_l2": float(np.linalg.norm(stress)),
            "stress_linf": float(np.max(abs_stress)),
            "stress_p50": _quantile(abs_stress, 0.50),
            "stress_p95": _quantile(abs_stress, 0.95),
            "stress_p99": _quantile(abs_stress, 0.99),
            "zero_mode_before_removal": full_field.last_zero_mode_magnitude,
            "mechanical_equilibrium_residual": full_field.last_equilibrium_residual,
            "source_increment_l2": float(np.linalg.norm(source)),
            "source_increment_linf": float(np.max(np.abs(source))),
            "source_increment_sum_00": float(np.sum(source[0, 0])),
            "source_increment_sum_01": float(np.sum(source[0, 1])),
            "source_increment_sum_10": float(np.sum(source[1, 0])),
            "source_increment_sum_11": float(np.sum(source[1, 1])),
            "interfacial_energy": interfacial,
            "elastic_energy": elastic, "total_energy": interfacial + elastic,
            "pf_interfacial_energy_change": interfacial - self.previous_interfacial,
            "source_elastic_energy_change": source_delta,
            "predicted_source_work": self.predicted_source_work,
            "source_work_error": source_work_error,
            "integrated_boundary_sweep": self.boundary_sweep,
            "integrated_source_measure": self.source_measure,
            "elastic_capillary_ratio_max": float(np.max(ratios)) if len(ratios) else 0.0,
            "elastic_capillary_ratio_p95": _quantile(ratios, 0.95) if len(ratios) else 0.0,
            "largest_one_step_population_loss": len(lost),
            "largest_100_step_population_loss": loss_100,
            "extinction_cluster_size": cluster,
            **metrics,
        }
        self.scalar_rows.append(row)
        self.previous_ids = current_ids
        self.previous_centroids = centroids
        self.previous_interfacial = interfacial
        self._evaluate_guard(row, simulation)
        if diag.step >= self.field_start and diag.step % self.field_cadence == 0:
            self.save_fields(simulation, "cadence")
        if self.guard_reason is not None:
            self.save_fields(simulation, "guard")
            self.flush()
        elif len(self.scalar_rows) >= self.flush_rows:
            self.flush()
        return self.guard_reason if self.terminate_on_guard else None

    def _evaluate_guard(self, row: dict[str, Any], simulation: Any) -> None:
        reasons: list[str] = []
        if row["largest_100_step_population_loss"] > 0.10 * max(
            row["grain_count"] + row["largest_100_step_population_loss"], 1
        ):
            reasons.append("grain_count_drop_gt_10pct_in_100_steps")
        if row["compactness_mean"] > 2.5:
            reasons.append("mean_compactness_gt_2.5")
        if row["compactness_max"] > 6.0:
            reasons.append("max_compactness_gt_6")
        if row["clipped_fraction"] > self.clip_spike_fraction:
            reasons.append("clipping_spike")
        if row["newly_extinct_phases"] >= self.extinction_spike:
            reasons.append("extinction_spike")
        for key in ("total_energy", "stress_linf", "eigenstrain_linf"):
            if not np.isfinite(row[key]):
                reasons.append(f"nonfinite_{key}")
        if not np.all(np.isfinite(simulation.driving_field)):
            reasons.append("nonfinite_driving")
        if reasons and self.guard_reason is None:
            self.guard_reason = ";".join(reasons)
            (self.output_dir / "diagnostic_capture.json").write_text(json.dumps({
                "schema_version": self.schema_version, "step": row["step"],
                "time": row["time"], "reasons": reasons,
                "scientific_status": "diagnostic_capture_not_completed",
            }, indent=2) + "\n")

    def save_fields(self, simulation: Any, reason: str) -> Path:
        path = self.field_dir / f"step-{simulation.solver.step_number:07d}-{reason}.npz"
        np.savez_compressed(
            path, eta=simulation.solver.eta,
            active_phases=simulation.solver.active_phases,
            labels=simulation.solver.labels,
            eigenstrain=simulation.full_field.eigenstrain,
            stress=simulation.full_field.stress,
            source_increment=simulation.full_field.source_increment,
            driving_field=simulation.driving_field,
            step=np.asarray(simulation.solver.step_number),
            time=np.asarray(simulation.solver.time),
        )
        return path

    def flush(self) -> None:
        if self.scalar_rows:
            table = pa.Table.from_pylist(self.scalar_rows)
            pq.write_table(table, self.scalar_dir / f"part-{self.scalar_part:06d}.parquet", compression="zstd")
            self.scalar_part += 1
            self.scalar_rows.clear()
        if self.boundary_rows:
            table = pa.Table.from_pylist(self.boundary_rows)
            pq.write_table(table, self.boundary_dir / f"part-{self.boundary_part:06d}.parquet", compression="zstd")
            self.boundary_part += 1
            self.boundary_rows.clear()

    def checkpoint_state(self) -> dict[str, Any]:
        self.flush()
        return {
            "schema_version": self.schema_version,
            "scalar_part": self.scalar_part, "boundary_part": self.boundary_part,
            "population_history": list(self.population_history),
            "previous_ids": sorted(self.previous_ids),
            "previous_centroids": {str(key): value for key, value in self.previous_centroids.items()},
            "previous_interfacial": self.previous_interfacial,
            "guard_reason": self.guard_reason,
        }

    def restore(self, state: dict[str, Any]) -> None:
        if int(state.get("schema_version", -1)) != self.schema_version:
            raise ValueError("unsupported QIU diagnostic checkpoint schema")
        self.scalar_part = int(state["scalar_part"])
        self.boundary_part = int(state["boundary_part"])
        self.population_history = deque((int(a), int(b)) for a, b in state["population_history"])
        self.previous_ids = {int(value) for value in state["previous_ids"]}
        self.previous_centroids = {
            int(key): tuple(map(float, value)) for key, value in state["previous_centroids"].items()
        }
        self.previous_interfacial = float(state["previous_interfacial"])
        self.guard_reason = state.get("guard_reason")

    def close(self) -> None:
        self.flush()
