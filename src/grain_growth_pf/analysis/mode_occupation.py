from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from grain_growth_pf.io.event_ledger import event_ledger_path, iter_event_ledger
from grain_growth_pf.io.provenance import git_sha


MODE_COMPLETION_TYPES = {"disconnection_mode", "compatibility_release"}


def _run_mode_occupation(run_dir: Path) -> dict[str, Any]:
    manifest = json.loads((run_dir / "manifest.json").read_text())
    count = 0
    signed_sum = 0.0
    absolute_sum = 0.0
    square_sum = 0.0
    shear_strain_sum = 0.0
    volumetric_strain_sum = 0.0
    family_counts: Counter[str] = Counter()
    magnitude_counts: Counter[str] = Counter()
    columns = [
        "event_type", "barrier_type", "activation_volume",
        "shear_strain_increment", "volumetric_strain_increment",
    ]
    for frame in iter_event_ledger(event_ledger_path(run_dir), columns=columns):
        completed = frame[frame["event_type"].isin(MODE_COMPLETION_TYPES)]
        if completed.empty:
            continue
        components = completed["activation_volume"].astype("string").str.split(
            ";", n=1, expand=True
        )
        if components.shape[1] != 2:
            continue
        normal = pd.to_numeric(components[0], errors="coerce").to_numpy(float)
        shear = pd.to_numeric(components[1], errors="coerce").to_numpy(float)
        valid = np.isfinite(normal) & np.isfinite(shear) & (normal != 0)
        if not np.any(valid):
            continue
        beta = shear[valid] / normal[valid]
        absolute_beta = np.abs(beta)
        count += int(len(beta))
        signed_sum += float(beta.sum())
        absolute_sum += float(absolute_beta.sum())
        square_sum += float(np.square(beta).sum())
        families = completed.loc[valid, "barrier_type"].dropna().astype(str)
        family_counts.update(families)
        magnitude_counts.update(f"{value:.12g}" for value in absolute_beta)
        shear_strain_sum += float(pd.to_numeric(
            completed.loc[valid, "shear_strain_increment"], errors="coerce"
        ).fillna(0.0).sum())
        volumetric_strain_sum += float(pd.to_numeric(
            completed.loc[valid, "volumetric_strain_increment"], errors="coerce"
        ).fillna(0.0).sum())
    return {
        "run": str(run_dir),
        "seed": int(manifest["config"]["seed"]),
        "temperature_K": float(manifest["config"]["pf"]["temperature"]),
        "simulation_git_sha": manifest["git_sha"],
        "mode_events": count,
        "signed_beta_sum": signed_sum,
        "absolute_beta_sum": absolute_sum,
        "squared_beta_sum": square_sum,
        "mean_signed_beta": signed_sum / count if count else None,
        "mean_absolute_beta": absolute_sum / count if count else None,
        "family_counts": dict(sorted(family_counts.items())),
        "absolute_beta_counts": dict(sorted(magnitude_counts.items(), key=lambda item: float(item[0]))),
        "signed_shear_strain": shear_strain_sum,
        "signed_volumetric_strain": volumetric_strain_sum,
    }


def _pooled_summary(runs: list[dict[str, Any]], regime: str, temperature: float,
                    bootstrap_samples: int) -> dict[str, Any]:
    count = sum(run["mode_events"] for run in runs)
    signed_sum = sum(run["signed_beta_sum"] for run in runs)
    absolute_sum = sum(run["absolute_beta_sum"] for run in runs)
    square_sum = sum(run["squared_beta_sum"] for run in runs)
    family_counts: Counter[str] = Counter()
    magnitude_counts: Counter[str] = Counter()
    for run in runs:
        family_counts.update(run["family_counts"])
        magnitude_counts.update(run["absolute_beta_counts"])
    digest = hashlib.sha256(f"mode-occupation:{regime}:{temperature}".encode()).digest()
    rng = np.random.default_rng(int.from_bytes(digest[:8], "little"))
    bootstrap_absolute, bootstrap_signed = [], []
    for _ in range(bootstrap_samples):
        selected = rng.integers(0, len(runs), len(runs))
        sampled = [runs[index] for index in selected]
        sampled_count = sum(run["mode_events"] for run in sampled)
        if sampled_count:
            bootstrap_absolute.append(
                sum(run["absolute_beta_sum"] for run in sampled) / sampled_count
            )
            bootstrap_signed.append(
                sum(run["signed_beta_sum"] for run in sampled) / sampled_count
            )
    variance = square_sum / count - (signed_sum / count) ** 2 if count else np.nan
    return {
        "temperature_K": temperature,
        "realizations": len(runs),
        "mode_events": count,
        "mean_signed_beta": signed_sum / count if count else None,
        "mean_signed_beta_95pct": (
            np.quantile(bootstrap_signed, [0.025, 0.975]).tolist()
            if bootstrap_signed else None
        ),
        "mean_absolute_beta": absolute_sum / count if count else None,
        "mean_absolute_beta_95pct": (
            np.quantile(bootstrap_absolute, [0.025, 0.975]).tolist()
            if bootstrap_absolute else None
        ),
        "beta_standard_deviation": float(np.sqrt(max(0.0, variance))) if count else None,
        "family_counts": dict(sorted(family_counts.items())),
        "family_fractions": {
            key: value / count for key, value in sorted(family_counts.items())
        } if count else {},
        "absolute_beta_counts": dict(sorted(
            magnitude_counts.items(), key=lambda item: float(item[0])
        )),
        "absolute_beta_fractions": {
            key: value / count for key, value in sorted(
                magnitude_counts.items(), key=lambda item: float(item[0])
            )
        } if count else {},
        "signed_shear_strain": sum(run["signed_shear_strain"] for run in runs),
        "signed_volumetric_strain": sum(
            run["signed_volumetric_strain"] for run in runs
        ),
        "bootstrap_samples": bootstrap_samples,
    }


def analyze_mode_occupation(campaign_dir: str | Path, regime: str,
                            bootstrap_samples: int = 5000) -> dict[str, Any]:
    """Measure completed-mode occupation and effective beta by temperature."""
    campaign_dir = Path(campaign_dir)
    campaign = json.loads((campaign_dir / "campaign_manifest.json").read_text())
    runs = []
    for raw_run in campaign["runs"]:
        run_dir = Path(raw_run)
        manifest_path = run_dir / "manifest.json"
        if not manifest_path.exists():
            continue
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("status") != "completed":
            continue
        if manifest["config"]["regime"] != regime:
            continue
        runs.append(_run_mode_occupation(run_dir))
    if not runs:
        raise ValueError(f"no completed runs for regime {regime!r} in {campaign_dir}")
    grouped: dict[float, list[dict[str, Any]]] = {}
    for run in runs:
        grouped.setdefault(run["temperature_K"], []).append(run)
    temperatures = [
        _pooled_summary(grouped[temperature], regime, temperature, bootstrap_samples)
        for temperature in sorted(grouped)
    ]
    absolute = [item["mean_absolute_beta"] for item in temperatures]
    signed = [item["mean_signed_beta"] for item in temperatures]
    return {
        "status": "completed",
        "campaign": str(campaign_dir),
        "regime": regime,
        "analysis_git_sha": git_sha(),
        "simulation_git_shas": sorted({run["simulation_git_sha"] for run in runs}),
        "mode_completion_event_types": sorted(MODE_COMPLETION_TYPES),
        "estimator_note": (
            "Event-conditioned occupation over completed mode rows. beta is derived "
            "from activation_volume_shear/activation_volume_normal; it is not a "
            "time-weighted mobility or causal resistance fraction."
        ),
        "temperature_summaries": temperatures,
        "temperature_trend": {
            "mean_absolute_beta_monotone_increasing": bool(
                all(right >= left for left, right in zip(absolute, absolute[1:]))
            ),
            "mean_signed_beta_monotone_increasing": bool(
                all(right >= left for left, right in zip(signed, signed[1:]))
            ),
            "absolute_beta_high_to_low_ratio": (
                absolute[-1] / absolute[0] if absolute[0] else None
            ),
        },
        "run_summaries": runs,
    }
