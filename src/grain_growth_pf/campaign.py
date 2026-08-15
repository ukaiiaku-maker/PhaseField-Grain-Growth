from __future__ import annotations

import json
import multiprocessing as mp
import time
import traceback
from copy import deepcopy
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

from grain_growth_pf.config import ModelConfig, PFConfig
from grain_growth_pf.io.provenance import canonical_hash, git_sha, write_manifest
from grain_growth_pf.pf.initial_conditions import initial_condition_identity, prepare_initial_condition
from grain_growth_pf.simulation import EventResolvedSimulation


def _run_one(payload: tuple[dict[str, Any], str, bool, str]) -> dict[str, str]:
    data, path, resume, code_sha = payload
    try:
        config = ModelConfig.from_dict(data)
        EventResolvedSimulation(config, path, resume=resume, code_sha=code_sha).run()
        return {"path": path, "status": "completed"}
    except Exception:
        Path(path).mkdir(parents=True, exist_ok=True)
        failure = traceback.format_exc()
        (Path(path) / "traceback.log").write_text(failure)
        manifest_path = Path(path) / "manifest.json"
        previous = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
        write_manifest(manifest_path, data, "failed", {"failure": failure.splitlines()[-1]},
                       code_sha=previous.get("git_sha", git_sha()))
        return {"path": path, "status": "failed"}


def _prepare_one(payload: tuple[PFConfig, int, dict[str, Any], str, str]) -> str:
    pf, seed, parameters, path, code_sha = payload
    return str(prepare_initial_condition(pf, seed, parameters, path, code_sha))


def _run_index(root: Path) -> tuple[dict[str, str], dict[str, str]]:
    completed: dict[str, str] = {}
    resumable: dict[str, str] = {}
    if not root.exists():
        return completed, resumable
    for path in root.rglob("manifest.json"):
        try:
            manifest = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if manifest.get("status") == "completed" and manifest.get("config_sha256"):
            completed[manifest["config_sha256"]] = str(path.parent)
        checkpoint = path.parent / "checkpoint.json"
        arrays = path.parent / "checkpoint.npz"
        stale = checkpoint.exists() and time.time() - checkpoint.stat().st_mtime > 120
        if (manifest.get("status") == "failed" or stale) and checkpoint.exists() and arrays.exists() \
                and manifest.get("config_sha256"):
            resumable[manifest["config_sha256"]] = str(path.parent)
    return completed, resumable


def enumerate_campaign(spec: dict[str, Any]) -> list[ModelConfig]:
    base = ModelConfig.from_dict(spec.get("base", {}))
    regimes = spec["regimes"]
    temperatures = spec.get("temperatures", [base.pf.temperature])
    seeds = spec.get("seeds", [base.seed])
    result = []
    for regime in regimes:
        overrides = regimes[regime] or {}
        compatibility = overrides.get("compatibility_model", base.compatibility_model)
        if compatibility is False:
            compatibility = "off"
        mechanics = overrides.get("mechanics_backend", base.mechanics_backend)
        modules = tuple(overrides.get("modules", base.active_modules))
        parameters = {**base.parameters, **overrides.get("parameters", {})}
        pf_overrides = dict(overrides.get("pf", {}))
        pf_overrides.pop("temperature", None)
        if "shape" in pf_overrides:
            pf_overrides["shape"] = tuple(pf_overrides["shape"])
        for temperature in temperatures:
            for seed in seeds:
                pf = replace(base.pf, **pf_overrides, temperature=float(temperature))
                result.append(replace(base, regime=regime, seed=int(seed), pf=pf,
                                      active_modules=modules, compatibility_model=compatibility,
                                      mechanics_backend=mechanics, parameters=parameters,
                                      output_cadence=int(overrides.get("output_cadence", base.output_cadence)),
                                      max_steps=int(overrides.get("max_steps", base.max_steps)),
                                      termination_grains=int(overrides.get(
                                          "termination_grains", base.termination_grains))))
    return result


def launch_campaign(spec_path: str | Path, root: str | Path = "results/campaigns",
                    processes: int | None = None) -> Path:
    with Path(spec_path).open(encoding="utf-8") as handle:
        spec = yaml.safe_load(handle)
    submitted_spec = deepcopy(spec)
    if "regime_catalog" in spec:
        catalog_path = (Path(spec_path).parent / spec["regime_catalog"]).resolve()
        with catalog_path.open(encoding="utf-8") as handle:
            catalog = yaml.safe_load(handle)["regimes"]
        names = spec.pop("regime_names", list(catalog))
        spec["regimes"] = {name: catalog[name] for name in names}
    configs = enumerate_campaign(spec)
    code_sha = git_sha()
    initial_condition_files: dict[int, str] = {}
    if spec.get("prepare_initial_conditions", False):
        base = ModelConfig.from_dict(spec.get("base", {}))
        cache_root = Path(root).parent / "initial_conditions"
        preparation_payloads = []
        for seed in sorted({config.seed for config in configs}):
            identity_hash = initial_condition_identity(base.pf, seed, base.parameters, code_sha)
            state_path = cache_root / f"seed-{seed}-{identity_hash[:16]}.npz"
            initial_condition_files[seed] = str(state_path)
            preparation_payloads.append((base.pf, seed, base.parameters, str(state_path), code_sha))
        preparation_workers = min(
            processes or max(1, mp.cpu_count() - 1), len(preparation_payloads)
        ) if preparation_payloads else 0
        if preparation_workers == 1:
            [_prepare_one(payload) for payload in preparation_payloads]
        elif preparation_workers > 1:
            with mp.get_context("spawn").Pool(preparation_workers) as pool:
                pool.map(_prepare_one, preparation_payloads)
        configs = [replace(
            config, parameters={**config.parameters, "initial_state_file": initial_condition_files[config.seed]}
        ) for config in configs]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    identity = canonical_hash({"runs": [c.to_dict() for c in configs]})[:10]
    campaign_dir = Path(root) / f"{stamp}-{identity}"
    campaign_dir.mkdir(parents=True, exist_ok=False)
    payloads = []
    reused: list[str] = []
    resumed: list[str] = []
    prior, resumable = _run_index(Path(root))
    for config in configs:
        full_hash = canonical_hash({"config": config.to_dict(), "git_sha": code_sha})
        run_hash = full_hash[:12]
        if full_hash in prior:
            reused.append(prior[full_hash])
            continue
        if full_hash in resumable:
            resumed.append(resumable[full_hash])
            payloads.append((config.to_dict(), resumable[full_hash], True, code_sha))
        else:
            run_dir = campaign_dir / f"{config.regime}-T{config.pf.temperature:g}-s{config.seed}-{run_hash}"
            payloads.append((config.to_dict(), str(run_dir), False, code_sha))
    (campaign_dir / "campaign_manifest.json").write_text(json.dumps({
        "source_spec": str(spec_path), "runs": reused + [p[1] for p in payloads],
        "specification": submitted_spec,
        "initial_condition_files": initial_condition_files,
        "reused_completed": reused, "resumed_runs": resumed, "status": "running"
    }, indent=2) + "\n")
    workers = min(processes or max(1, mp.cpu_count() - 1), len(payloads)) if payloads else 0
    if workers == 0:
        outcomes = []
    elif workers == 1:
        outcomes = [_run_one(p) for p in payloads]
    else:
        with mp.get_context("spawn").Pool(workers) as pool:
            outcomes = pool.map(_run_one, payloads)
    completed = reused + [o["path"] for o in outcomes]
    failures = [o["path"] for o in outcomes if o["status"] != "completed"]
    (campaign_dir / "campaign_manifest.json").write_text(json.dumps({
        "source_spec": str(spec_path), "runs": completed,
        "specification": submitted_spec,
        "initial_condition_files": initial_condition_files,
        "status": "failed" if failures else "completed", "workers": workers,
        "reused_completed": reused, "resumed_runs": resumed, "failed_runs": failures,
    }, indent=2) + "\n")
    return campaign_dir
