from __future__ import annotations

import json
import multiprocessing as mp
import traceback
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

from grain_growth_pf.config import ModelConfig, PFConfig
from grain_growth_pf.io.provenance import canonical_hash, git_sha
from grain_growth_pf.simulation import EventResolvedSimulation


def _run_one(payload: tuple[dict[str, Any], str, bool]) -> dict[str, str]:
    data, path, resume = payload
    try:
        config = ModelConfig.from_dict(data)
        EventResolvedSimulation(config, path, resume=resume).run()
        return {"path": path, "status": "completed"}
    except Exception:
        Path(path).mkdir(parents=True, exist_ok=True)
        (Path(path) / "traceback.log").write_text(traceback.format_exc())
        return {"path": path, "status": "failed"}


def _completed_index(root: Path) -> dict[str, str]:
    completed: dict[str, str] = {}
    if not root.exists():
        return completed
    for path in root.rglob("manifest.json"):
        try:
            manifest = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if manifest.get("status") == "completed" and manifest.get("config_sha256"):
            completed[manifest["config_sha256"]] = str(path.parent)
    return completed


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
        for temperature in temperatures:
            for seed in seeds:
                pf = replace(base.pf, temperature=float(temperature))
                result.append(replace(base, regime=regime, seed=int(seed), pf=pf,
                                      active_modules=modules, compatibility_model=compatibility,
                                      mechanics_backend=mechanics))
    return result


def launch_campaign(spec_path: str | Path, root: str | Path = "results/campaigns",
                    processes: int | None = None) -> Path:
    with Path(spec_path).open(encoding="utf-8") as handle:
        spec = yaml.safe_load(handle)
    if "regime_catalog" in spec:
        catalog_path = (Path(spec_path).parent / spec["regime_catalog"]).resolve()
        with catalog_path.open(encoding="utf-8") as handle:
            catalog = yaml.safe_load(handle)["regimes"]
        names = spec.pop("regime_names", list(catalog))
        spec["regimes"] = {name: catalog[name] for name in names}
    configs = enumerate_campaign(spec)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    identity = canonical_hash({"runs": [c.to_dict() for c in configs]})[:10]
    campaign_dir = Path(root) / f"{stamp}-{identity}"
    campaign_dir.mkdir(parents=True, exist_ok=False)
    payloads = []
    reused: list[str] = []
    prior = _completed_index(Path(root))
    code_sha = git_sha()
    for config in configs:
        full_hash = canonical_hash({"config": config.to_dict(), "git_sha": code_sha})
        run_hash = full_hash[:12]
        if full_hash in prior:
            reused.append(prior[full_hash])
            continue
        run_dir = campaign_dir / f"{config.regime}-T{config.pf.temperature:g}-s{config.seed}-{run_hash}"
        payloads.append((config.to_dict(), str(run_dir), False))
    (campaign_dir / "campaign_manifest.json").write_text(json.dumps({
        "source_spec": str(spec_path), "runs": reused + [p[1] for p in payloads],
        "reused_completed": reused, "status": "running"
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
        "status": "failed" if failures else "completed", "workers": workers,
        "reused_completed": reused, "failed_runs": failures,
    }, indent=2) + "\n")
    return campaign_dir
