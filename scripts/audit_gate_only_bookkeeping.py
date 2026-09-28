#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from grain_growth_pf.io.event_ledger import event_ledger_path, iter_event_ledger


def completed_runs(campaign: Path) -> list[Path]:
    manifest = json.loads((campaign / "campaign_manifest.json").read_text())
    runs = []
    for raw in manifest.get("runs", []):
        run = Path(raw)
        mpath = run / "manifest.json"
        if not mpath.exists():
            continue
        data = json.loads(mpath.read_text())
        if data.get("status") == "completed":
            runs.append(run)
    return runs


def audit_run(run: Path) -> dict[str, object]:
    manifest = json.loads((run / "manifest.json").read_text())
    config = manifest["config"]
    closure = str(config.get("parameters", {}).get("migration_closure", ""))
    checkpoint = json.loads((run / "checkpoint.json").read_text())
    ledgers = []
    free_volume_violations = 0
    for state in checkpoint.get("domains", {}).values():
        ledgers.append(abs(float(state.get("normal_displacement_ledger", 0.0))))
        ledgers.append(abs(float(state.get("normal_release_remaining", 0.0))))
        fv = state.get("free_volume", {})
        required = float(fv.get("required_total", 0.0))
        accommodated = float(fv.get("accommodated_total", 0.0))
        if required < -1e-12 or accommodated < -1e-12 or accommodated > required + 1e-10:
            free_volume_violations += 1

    climb_rows = 0
    climb_release_q = 0.0
    climb_nv = 0.0
    path = event_ledger_path(run)
    if path.exists():
        for frame in iter_event_ledger(
            path, columns=["event_type", "release_Delta_q", "Nv"]
        ):
            if "event_type" not in frame:
                continue
            mask = frame["event_type"] == "climb_quota_completion"
            if not mask.any():
                continue
            rows = frame.loc[mask]
            climb_rows += len(rows)
            climb_release_q += float(pd.to_numeric(
                rows["release_Delta_q"], errors="coerce"
            ).fillna(0.0).abs().sum())
            climb_nv += float(pd.to_numeric(
                rows["Nv"], errors="coerce"
            ).fillna(0.0).abs().sum())

    max_hidden = max(ledgers) if ledgers else 0.0
    passed = True
    reasons = []
    if closure == "gate_only" and max_hidden > 1e-10:
        passed = False
        reasons.append(f"gate-only hidden PF displacement remains: {max_hidden:g}")
    if closure == "gate_only" and climb_rows and (climb_release_q > 1e-10 or climb_nv > 1e-10):
        passed = False
        reasons.append(
            f"climb summary applied second q increment: release={climb_release_q:g}, Nv={climb_nv:g}"
        )
    if free_volume_violations:
        passed = False
        reasons.append(f"free-volume balance violations={free_volume_violations}")
    return {
        "run": str(run),
        "regime": config["regime"],
        "migration_closure": closure,
        "max_hidden_normal_pf_state": max_hidden,
        "climb_completion_rows": climb_rows,
        "climb_completion_abs_release_Delta_q": climb_release_q,
        "climb_completion_abs_Nv": climb_nv,
        "free_volume_violations": free_volume_violations,
        "PASS": passed,
        "reason": "; ".join(reasons),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    rows = [audit_run(run) for run in completed_runs(Path(args.campaign))]
    if not rows:
        raise SystemExit("no completed runs to audit")
    table = pd.DataFrame(rows)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, index=False)
    print(table.to_string(index=False))
    if not bool(table["PASS"].all()):
        raise SystemExit("gate-only bookkeeping audit FAILED")
    print("gate-only bookkeeping audit PASS")


if __name__ == "__main__":
    main()
