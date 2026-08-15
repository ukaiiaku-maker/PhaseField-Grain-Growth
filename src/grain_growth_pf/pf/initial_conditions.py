from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from grain_growth_pf.config import PFConfig
from grain_growth_pf.io.provenance import canonical_hash
from grain_growth_pf.pf.geometry import voronoi_polycrystal
from grain_growth_pf.pf.solver import MultiphaseFieldSolver


def initial_condition_identity(pf: PFConfig, seed: int, parameters: dict[str, Any],
                               code_sha: str) -> str:
    controls = {
        "seed": seed, "shape": pf.shape, "grid_spacing": pf.grid_spacing,
        "interface_width": pf.interface_width, "time_step": pf.time_step,
        "gb_energy": pf.gb_energy, "intrinsic_mobility": pf.intrinsic_mobility,
        "boundary_conditions": pf.boundary_conditions,
        "grain_extinction_threshold": pf.grain_extinction_threshold,
        "initial_grains": parameters.get("initial_grains", 50),
        "equilibration_steps": parameters.get("equilibration_steps", 0),
        "equilibrate_to_grains": parameters.get("equilibrate_to_grains"),
        "equilibration_max_steps": parameters.get("equilibration_max_steps", 5000),
        "code_sha": code_sha,
    }
    return canonical_hash(controls)


def prepare_initial_condition(pf: PFConfig, seed: int, parameters: dict[str, Any],
                              path: str | Path, code_sha: str) -> Path:
    """Create a compact, pre-coarsened state shared across mechanisms and temperatures."""
    target = Path(path)
    metadata_path = target.with_suffix(".json")
    if target.exists() and metadata_path.exists():
        metadata = json.loads(metadata_path.read_text())
        if metadata.get("status") == "completed":
            return target
    target.parent.mkdir(parents=True, exist_ok=True)
    eta, seed_positions, orientations = voronoi_polycrystal(
        pf.shape, int(parameters.get("initial_grains", 50)), seed,
        width=pf.interface_width / 2, periodic=pf.boundary_conditions == "periodic",
    )
    solver = MultiphaseFieldSolver(eta, pf)
    steps = int(parameters.get("equilibration_steps", 0))
    for _ in range(steps):
        solver.step()
    target_grains = parameters.get("equilibrate_to_grains")
    maximum = int(parameters.get("equilibration_max_steps", 5000))
    if target_grains is not None:
        desired = int(target_grains)
        while np.count_nonzero(solver.active_phases) > desired and steps < maximum:
            solver.step()
            steps += 1
            if steps % 100 == 0:
                metadata_path.write_text(json.dumps({
                    "status": "running", "steps": steps,
                    "active_grains": int(np.count_nonzero(solver.active_phases)),
                    "target_grains": desired, "git_sha": code_sha,
                }, indent=2) + "\n")
        if np.count_nonzero(solver.active_phases) > desired:
            metadata_path.write_text(json.dumps({
                "status": "failed", "steps": steps,
                "active_grains": int(np.count_nonzero(solver.active_phases)),
                "target_grains": desired, "git_sha": code_sha,
            }, indent=2) + "\n")
            raise RuntimeError(
                f"cached pre-equilibration retained {np.count_nonzero(solver.active_phases)} "
                f"grains after {maximum} steps"
            )
    active_original_ids = np.flatnonzero(solver.active_phases)
    eta = solver.eta[active_original_ids].copy()
    orientations = orientations[active_original_ids].copy()
    np.savez_compressed(
        target, eta=eta, orientations=orientations, seed_positions=seed_positions,
        active_original_ids=active_original_ids,
        equilibration_steps=np.asarray(steps),
    )
    metadata_path.write_text(json.dumps({
        "status": "completed", "steps": steps, "active_grains": len(active_original_ids),
        "target_grains": target_grains, "git_sha": code_sha,
        "identity": initial_condition_identity(pf, seed, parameters, code_sha),
    }, indent=2) + "\n")
    return target
