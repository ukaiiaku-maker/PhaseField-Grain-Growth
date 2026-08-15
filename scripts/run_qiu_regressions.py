#!/usr/bin/env python3
from __future__ import annotations

import json
import platform
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

from grain_growth_pf.config import ModelConfig, PFConfig
from grain_growth_pf.io.provenance import git_sha
from grain_growth_pf.simulation import EventResolvedSimulation


def _boundary_statistics(path: Path) -> dict[str, float]:
    tracks = pd.read_csv(path / "boundary_tracks.csv")
    curvature = tracks["curvature"].to_numpy(float)
    velocity = tracks["normal_velocity"].to_numpy(float)
    valid = np.isfinite(curvature) & np.isfinite(velocity)
    valid &= (np.abs(curvature) > 1e-12) & (np.abs(velocity) > 1e-12)
    if np.count_nonzero(valid) > 2:
        correlation = float(np.corrcoef(curvature[valid], velocity[valid])[0, 1])
        reverse = float(np.mean(curvature[valid] * velocity[valid] < 0.0))
    else:
        correlation, reverse = np.nan, np.nan
    return {
        "velocity_curvature_correlation": correlation,
        "reverse_curvature_fraction": reverse,
        "boundary_samples": int(np.count_nonzero(valid)),
    }


def _matched_pair(name: str, shape: tuple[int, int], grains: int, steps: int,
                  seed: int, root: Path, code_sha: str) -> dict:
    common = ModelConfig(
        regime=f"{name}-curvature",
        seed=seed,
        pf=PFConfig(
            shape=shape, interface_width=4.0, time_step=0.03,
            intrinsic_mobility=0.5, adaptive_stepping=True,
        ),
        output_cadence=2,
        max_steps=steps,
        termination_grains=1,
        parameters={"initial_grains": grains, "easy_beta": 0.8},
    )
    shear = replace(
        common,
        regime=f"{name}-qiu-shear",
        mechanics_backend="qiu_full_field",
        active_modules=("qiu_reference_shear",),
    )
    control_simulation = EventResolvedSimulation(common, root / f"{name}-curvature", code_sha=code_sha)
    control_simulation.run()
    shear_simulation = EventResolvedSimulation(shear, root / f"{name}-qiu-shear", code_sha=code_sha)
    shear_simulation.run()
    assert shear_simulation.full_field is not None
    label_difference = float(np.mean(
        control_simulation.solver.labels != shear_simulation.solver.labels
    ))
    result = {
        "geometry": name,
        "shape": list(shape),
        "initial_grains": grains,
        "steps": steps,
        "control": _boundary_statistics(root / f"{name}-curvature"),
        "qiu_shear": _boundary_statistics(root / f"{name}-qiu-shear"),
        "eigenstrain_l2": float(np.linalg.norm(shear_simulation.full_field.eigenstrain)),
        "stress_l2": float(np.linalg.norm(shear_simulation.full_field.stress)),
        "feedback_l2": float(np.linalg.norm(shear_simulation.driving_field)),
        "label_difference_fraction": label_difference,
        "accumulated_shear_strain": shear_simulation.accumulated_shear_strain,
    }
    result["validation_passed"] = bool(
        result["eigenstrain_l2"] > 0.0
        and result["stress_l2"] > 0.0
        and result["feedback_l2"] > 0.0
        and result["label_difference_fraction"] > 0.0
        and np.isfinite(list(result["qiu_shear"].values())).all()
    )
    return result


def main() -> None:
    code_sha = git_sha()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    raw_root = Path("results/runs") / f"qiu-regressions-{stamp}-{code_sha[:8]}"
    raw_root.mkdir(parents=True, exist_ok=False)
    intrinsic = json.loads(Path("results/validation/numerical_validation.json").read_text())
    results = [
        _matched_pair("four-grain", (48, 48), 4, 160, 8101, raw_root, code_sha),
        _matched_pair("polycrystal", (64, 64), 24, 220, 8102, raw_root, code_sha),
    ]
    report = {
        "git_sha": code_sha,
        "raw_run_directory": str(raw_root),
        "intrinsic_curvature_regression": {
            "source": "results/validation/numerical_validation.json",
            "source_git_sha": intrinsic["git_sha"],
            "passed": intrinsic["passed"],
            "circle_max_relative_error": max(
                item["relative_error"] for item in intrinsic["circular_grains"]
            ),
        },
        "shear_regressions": results,
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
    }
    report["validation_passed"] = bool(
        report["intrinsic_curvature_regression"]["passed"]
        and all(item["validation_passed"] for item in results)
    )
    target = Path("results/validation/qiu_regression_benchmarks.json")
    target.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
