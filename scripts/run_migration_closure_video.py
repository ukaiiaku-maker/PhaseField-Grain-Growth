#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from grain_growth_pf.campaign import enumerate_campaign
from grain_growth_pf.config import ModelConfig
from grain_growth_pf.io.provenance import canonical_hash, git_sha
from grain_growth_pf.migration_closure import MigrationClosureSimulation
from grain_growth_pf.pf.initial_conditions import initial_condition_identity, prepare_initial_condition


class ClosureFrameSimulation(MigrationClosureSimulation):
    def _grain_size_moments(self) -> tuple[float, float, float]:
        areas = np.asarray(
            [grain.area for grain in self.snapshot.grains.values()], dtype=float
        )
        if not len(areas) or float(areas.sum()) <= 0.0:
            return float("nan"), float("nan"), float("nan")
        diameters = 2.0 * np.sqrt(areas / np.pi)
        mean = float(diameters.mean())
        population = float(2.0 * np.sqrt(areas.sum() / (np.pi * len(areas))))
        area_weighted = float(np.average(diameters, weights=areas))
        return mean, population, area_weighted

    def _write_frame(self, force: bool = False) -> None:
        cadence = max(1, int(self.config.parameters.get("video_frame_cadence", 10)))
        step = int(self.solver.step_number)
        frame_dir = self.output_dir / "frames"
        frame_dir.mkdir(parents=True, exist_ok=True)
        path = frame_dir / f"frame-{step:07d}.npz"
        if path.exists():
            return
        g_mean, g_population, g_area_weighted = self._grain_size_moments()
        progress_fraction = float(
            self.config.parameters.get("video_frame_progress_fraction", 0.0)
        )
        if not hasattr(self, "_initial_video_characteristic_size"):
            existing = sorted(frame_dir.glob("frame-*.npz"))
            if existing:
                with np.load(existing[0]) as first:
                    initial = float(first.get("G_population", g_population))
                with np.load(existing[-1]) as last:
                    previous = float(last.get("G_population", g_population))
            else:
                initial = previous = g_population
            self._initial_video_characteristic_size = initial
            self._last_video_characteristic_size = previous
        step_due = step % cadence == 0
        progress_due = bool(
            progress_fraction > 0.0
            and np.isfinite(g_population)
            and np.isfinite(self._last_video_characteristic_size)
            and g_population >= self._last_video_characteristic_size * (1.0 + progress_fraction)
        )
        if not force and not step_due and not progress_due:
            return
        labels = self.solver.labels.astype(np.uint16, copy=True)
        blocked = np.zeros(self.config.pf.shape, dtype=np.uint8)
        boundary_mask = np.zeros(self.config.pf.shape, dtype=np.uint8)
        shear = np.zeros(self.config.pf.shape, dtype=np.float32)
        shear_stress = np.zeros(self.config.pf.shape, dtype=np.float32)
        p_cap = np.zeros(self.config.pf.shape, dtype=np.float32)
        p_chem = np.zeros(self.config.pf.shape, dtype=np.float32)
        p_shear = np.zeros(self.config.pf.shape, dtype=np.float32)
        p_event = np.zeros(self.config.pf.shape, dtype=np.float32)
        p_net = np.zeros(self.config.pf.shape, dtype=np.float32)
        p_applied_total = np.zeros(self.config.pf.shape, dtype=np.float32)
        chi_s = np.full(self.config.pf.shape, np.nan, dtype=np.float32)
        local_shear_energy = np.zeros(self.config.pf.shape, dtype=np.float32)
        free_volume = np.zeros(self.config.pf.shape, dtype=np.float32)
        pending_state = np.zeros(self.config.pf.shape, dtype=np.uint8)
        climb_stage = np.zeros(self.config.pf.shape, dtype=np.uint8)
        sink_activity = np.zeros(self.config.pf.shape, dtype=np.uint8)
        stage_codes = {
            "inactive": 0, "nucleation": 1, "exchange": 2,
            "transport": 3, "quota_completion": 4,
        }
        shape = np.asarray(self.config.pf.shape)
        for segment in self.snapshot.boundaries.values():
            domain = self.domains.get(segment.entity_id)
            if domain is None or not len(segment.points):
                continue
            points = segment.points.astype(int) % shape
            yy, xx = points[:, 0], points[:, 1]
            boundary_mask[yy, xx] = 1
            if domain.blocked:
                blocked[yy, xx] = 1
            if domain.compatibility_pending:
                pending_state[yy, xx] |= 1
            if domain.area_loss_pending or domain.free_volume.deficit > 0.0:
                pending_state[yy, xx] |= 4
            shear[yy, xx] = float(domain.shear.state)
            shear_stress[yy, xx] = float(domain.shear.internal_shear_stress)
            balance = self._boundary_force_balance(domain, segment)
            p_cap[yy, xx] = balance.p_cap
            p_chem[yy, xx] = balance.p_chem
            p_shear[yy, xx] = balance.p_shear
            p_event[yy, xx] = balance.p_event
            p_net[yy, xx] = balance.p_net
            p_applied_total[yy, xx] = balance.p_applied_total
            chi_s[yy, xx] = balance.chi_s
            local_shear_energy[yy, xx] = float(domain.shear.energy)
            local_inventory = (
                self.defect_inventory.active_vacancy.get(domain.entity_id, 0.0)
                + self.defect_inventory.active_interstitial.get(domain.entity_id, 0.0)
                if self.area_loss_enabled else domain.free_volume.deficit
            )
            free_volume[yy, xx] = float(local_inventory)
            climb_stage[yy, xx] = stage_codes.get(domain.climb.stage.value, 0)
            if domain.entity_id in self._last_gb_sink_entities:
                sink_activity[yy, xx] = 1
        for key, tj in self.snapshot.triple_junctions.items():
            domain = self.tj_domains.get(key)
            if domain is None:
                continue
            y, x = np.rint(tj.position).astype(int) % shape
            if domain.blocked:
                pending_state[y, x] |= 2
            if domain.area_loss_pending:
                pending_state[y, x] |= 4
            climb_stage[y, x] = stage_codes.get(domain.climb.stage.value, 0)
            if key in self._last_tj_sink_entities:
                sink_activity[y, x] = 2
        qiu_shear_stress = (
            self.full_field.stress[0, 1].astype(np.float32, copy=True)
            if self.full_field is not None
            else np.zeros(self.config.pf.shape, dtype=np.float32)
        )
        on_boundary = boundary_mask.astype(bool)
        boundary_shear = shear[on_boundary].astype(float)
        nonzero_boundary = np.abs(boundary_shear) > 1e-12
        stored_shear_energy = float(sum(
            domain.shear.energy for domain in self.domains.values()
        ))
        domain_rows = [
            (segment, self.domains[segment.entity_id], self._boundary_force_balance(
                self.domains[segment.entity_id], segment
            ))
            for segment in self.snapshot.boundaries.values()
            if segment.entity_id in self.domains and segment.length > 0.0
        ]
        total_gb_length = float(sum(segment.length for segment, _, _ in domain_rows))
        active_rows = [
            row for row in domain_rows if abs(row[1].shear.state) > 1e-12
        ]
        active_length = float(sum(segment.length for segment, _, _ in active_rows))
        active_energies = np.asarray(
            [domain.shear.energy for _, domain, _ in active_rows], dtype=float
        )
        active_weights = np.asarray(
            [segment.length for segment, _, _ in active_rows], dtype=float
        )
        active_tau = np.asarray(
            [abs(domain.shear.internal_shear_stress) for _, domain, _ in active_rows],
            dtype=float,
        )
        active_p_shear = np.asarray(
            [abs(balance.p_shear) for _, _, balance in active_rows], dtype=float
        )
        active_chi = np.asarray(
            [balance.chi_s for _, _, balance in active_rows], dtype=float
        )
        finite_chi = np.isfinite(active_chi)

        def domain_quantile(values: np.ndarray, q: float) -> float:
            return float(np.quantile(values, q)) if values.size else float("nan")

        def length_mean(values: np.ndarray, weights: np.ndarray) -> float:
            return (
                float(np.average(values, weights=weights))
                if values.size and float(weights.sum()) > 0.0 else float("nan")
            )
        inventory = (
            self.defect_inventory.diagnostics()
            if self.area_loss_enabled else {
                "N_required": float(sum(d.free_volume.required_total for d in self.domains.values())),
                "N_accommodated_GB": float(sum(d.free_volume.accommodated_total for d in self.domains.values())),
                "N_accommodated_TJ": 0.0,
                "N_active_deficit": float(sum(d.free_volume.deficit for d in self.domains.values())),
                "conservation_residual": 0.0,
            }
        )
        gb_domains = [
            self.domains[key] for key in self.snapshot.boundaries if key in self.domains
        ]
        tj_domains = [
            self.tj_domains[key]
            for key in self.snapshot.triple_junctions if key in self.tj_domains
        ]
        all_domains = gb_domains + tj_domains
        g_occupancy = float(np.mean([domain.compatibility_pending for domain in gb_domains])) if gb_domains else 0.0
        t_occupancy = float(np.mean([domain.blocked and not domain.area_loss_pending for domain in tj_domains])) if tj_domains else 0.0
        c_occupancy = float(np.mean([domain.area_loss_pending for domain in all_domains])) if all_domains else 0.0
        sink_total = len(self._last_gb_sink_entities) + len(self._last_tj_sink_entities)
        np.savez_compressed(
            path,
            labels=labels,
            blocked=blocked,
            boundary_mask=boundary_mask,
            shear=shear,
            shear_stress=shear_stress,
            p_cap=p_cap,
            p_chem=p_chem,
            p_shear=p_shear,
            p_event=p_event,
            p_net=p_net,
            p_applied_total=p_applied_total,
            chi_s=chi_s,
            local_shear_energy=local_shear_energy,
            qiu_shear_stress=qiu_shear_stress,
            free_volume=free_volume,
            pending_state=pending_state,
            climb_stage=climb_stage,
            sink_activity=sink_activity,
            mobility=self.solver.mobility_scale.astype(np.float32),
            time=np.asarray(float(self.solver.time)),
            step=np.asarray(step),
            temperature=np.asarray(float(self.config.pf.temperature)),
            seed=np.asarray(int(self.config.seed)),
            shear_stiffness=np.asarray(
                float(self.config.parameters.get("shear_stiffness", 0.0))
            ),
            grain_count=np.asarray(len(self.snapshot.grains)),
            G_mean=np.asarray(g_mean),
            G_population=np.asarray(g_population),
            G_area_weighted=np.asarray(g_area_weighted),
            G_over_G0=np.asarray(
                g_population / self._initial_video_characteristic_size
                if self._initial_video_characteristic_size > 0.0 else np.nan
            ),
            G_occupancy=np.asarray(g_occupancy),
            T_occupancy=np.asarray(t_occupancy),
            C_occupancy=np.asarray(c_occupancy),
            GB_sink_fraction=np.asarray(
                len(self._last_gb_sink_entities) / sink_total if sink_total else 0.0
            ),
            TJ_sink_fraction=np.asarray(
                len(self._last_tj_sink_entities) / sink_total if sink_total else 0.0
            ),
            N_required=np.asarray(float(inventory["N_required"])),
            N_accommodated_GB=np.asarray(float(inventory["N_accommodated_GB"])),
            N_accommodated_TJ=np.asarray(float(inventory["N_accommodated_TJ"])),
            N_active_deficit=np.asarray(float(inventory["N_active_deficit"])),
            conservation_residual=np.asarray(float(inventory["conservation_residual"])),
            shear_state_max_abs=np.asarray(float(np.max(np.abs(shear)))),
            shear_state_rms=np.asarray(float(np.sqrt(np.mean(shear.astype(float) ** 2)))),
            boundary_shear_state_rms=np.asarray(
                float(np.sqrt(np.mean(boundary_shear**2))) if boundary_shear.size else 0.0
            ),
            nonzero_gb_length_fraction=np.asarray(
                float(np.mean(nonzero_boundary)) if boundary_shear.size else 0.0
            ),
            stored_shear_energy=np.asarray(stored_shear_energy),
            total_gb_length=np.asarray(total_gb_length),
            active_shear_domain_count=np.asarray(len(active_rows)),
            active_shear_length=np.asarray(active_length),
            active_shear_bearing_gb_fraction=np.asarray(
                active_length / total_gb_length if total_gb_length else 0.0
            ),
            stored_shear_energy_per_active_gb_length=np.asarray(
                stored_shear_energy / active_length if active_length else 0.0
            ),
            stored_shear_energy_per_active_domain=np.asarray(
                float(active_energies.mean()) if active_energies.size else 0.0
            ),
            local_shear_energy_p50=np.asarray(domain_quantile(active_energies, 0.50)),
            local_shear_energy_p90=np.asarray(domain_quantile(active_energies, 0.90)),
            local_shear_energy_p95=np.asarray(domain_quantile(active_energies, 0.95)),
            local_shear_energy_p99=np.asarray(domain_quantile(active_energies, 0.99)),
            mean_abs_tau_int_active_domain=np.asarray(
                float(active_tau.mean()) if active_tau.size else 0.0
            ),
            mean_abs_tau_int_active_length=np.asarray(
                length_mean(active_tau, active_weights)
            ),
            mean_abs_p_shear_active_domain=np.asarray(
                float(active_p_shear.mean()) if active_p_shear.size else 0.0
            ),
            mean_abs_p_shear_active_length=np.asarray(
                length_mean(active_p_shear, active_weights)
            ),
            chi_s_p50_active_domain=np.asarray(domain_quantile(active_chi[finite_chi], 0.50)),
            chi_s_p90_active_domain=np.asarray(domain_quantile(active_chi[finite_chi], 0.90)),
            chi_s_p95_active_domain=np.asarray(domain_quantile(active_chi[finite_chi], 0.95)),
            chi_s_p99_active_domain=np.asarray(domain_quantile(active_chi[finite_chi], 0.99)),
            mean_chi_s_active_domain=np.asarray(
                float(active_chi[finite_chi].mean()) if np.any(finite_chi) else float("nan")
            ),
            mean_chi_s_active_length=np.asarray(
                length_mean(active_chi[finite_chi], active_weights[finite_chi])
            ),
            fraction_active_gb_length_chi_s_near_one=np.asarray(
                float(active_weights[finite_chi & (active_chi > 0.8) & (active_chi < 1.2)].sum()
                      / active_weights[finite_chi].sum())
                if np.any(finite_chi) and active_weights[finite_chi].sum() else 0.0
            ),
            fraction_active_gb_length_chi_s_gt_0p8=np.asarray(
                float(active_weights[finite_chi & (active_chi > 0.8)].sum()
                      / active_weights[finite_chi].sum())
                if np.any(finite_chi) and active_weights[finite_chi].sum() else 0.0
            ),
            fraction_active_gb_length_chi_s_gt_one=np.asarray(
                float(active_weights[finite_chi & (active_chi > 1.0)].sum()
                      / active_weights[finite_chi].sum())
                if np.any(finite_chi) and active_weights[finite_chi].sum() else 0.0
            ),
            shear_stress_max_abs=np.asarray(float(np.max(np.abs(shear_stress)))),
            qiu_shear_stress_max_abs=np.asarray(float(np.max(np.abs(qiu_shear_stress)))),
        )
        self._last_video_characteristic_size = g_population

    def _save_checkpoint(self) -> None:
        super()._save_checkpoint()
        self._write_frame(force=False)

    def _write_tracks(self) -> None:
        super()._write_tracks()
        self._write_frame(force=False)

    def run(self) -> Path:
        self._write_frame(force=True)
        result = super().run()
        self._write_frame(force=True)
        return result


def _worker(payload: tuple[dict[str, Any], str, str]) -> dict[str, str]:
    config_dict, output, sha = payload
    try:
        ClosureFrameSimulation(ModelConfig.from_dict(config_dict), output, code_sha=sha).run()
        return {"path": output, "status": "completed"}
    except BaseException as exc:
        Path(output).mkdir(parents=True, exist_ok=True)
        (Path(output) / "video_failure.txt").write_text(
            f"{type(exc).__name__}: {exc}\n", encoding="utf-8"
        )
        return {"path": output, "status": "failed", "error": str(exc)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "spec", nargs="?", default="configs/production/migration_closure_video_c5.yaml"
    )
    parser.add_argument("--output-root", default="results/migration_closure_video")
    parser.add_argument("--processes", type=int, default=4)
    args = parser.parse_args()

    spec_path = Path(args.spec)
    spec = yaml.safe_load(spec_path.read_text(encoding="utf-8"))
    configs = enumerate_campaign(spec)
    sha = git_sha()
    output_root = Path(args.output_root)

    if spec.get("prepare_initial_conditions", False):
        base = ModelConfig.from_dict(spec.get("base", {}))
        cache_root = output_root.parent / "initial_conditions"
        files: dict[int, str] = {}
        for seed in sorted({config.seed for config in configs}):
            identity = initial_condition_identity(base.pf, seed, base.parameters, sha)
            state = cache_root / f"seed-{seed}-{identity[:16]}.npz"
            prepare_initial_condition(base.pf, seed, base.parameters, state, sha)
            files[seed] = str(state)
        configs = [replace(
            config, parameters={**config.parameters, "initial_state_file": files[config.seed]},
        ) for config in configs]

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    identity = canonical_hash({"configs": [c.to_dict() for c in configs], "sha": sha})[:10]
    root = output_root / f"{stamp}-{identity}"
    root.mkdir(parents=True, exist_ok=False)
    payloads = [
        (config.to_dict(), str(root / f"{config.regime}-T{config.pf.temperature:g}-s{config.seed}"), sha)
        for config in configs
    ]
    workers = min(max(1, args.processes), len(payloads))
    root_manifest = root / "video_manifest.json"
    root_manifest.write_text(json.dumps({
        "status": "running",
        "git_sha": sha,
        "source_spec": args.spec,
        "workers": workers,
        "runs": [{"path": payload[1], "status": "pending"} for payload in payloads],
    }, indent=2) + "\n", encoding="utf-8")
    if workers == 1:
        outcomes = [_worker(payload) for payload in payloads]
    else:
        with mp.get_context("spawn").Pool(workers) as pool:
            outcomes = pool.map(_worker, payloads, chunksize=1)
    status = "completed" if all(item["status"] == "completed" for item in outcomes) else "failed"
    root_manifest.write_text(json.dumps({
        "status": status,
        "git_sha": sha,
        "source_spec": args.spec,
        "workers": workers,
        "runs": outcomes,
    }, indent=2) + "\n", encoding="utf-8")
    print(root)
    if status != "completed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
