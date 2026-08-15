#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from pathlib import Path


def campaign_status(campaign: Path) -> dict[str, object]:
    manifest = json.loads((campaign / "campaign_manifest.json").read_text())
    rows = []
    for raw_run in manifest.get("runs", []):
        run = Path(raw_run)
        run_manifest_path = run / "manifest.json"
        checkpoint_path = run / "checkpoint.json"
        run_manifest = (
            json.loads(run_manifest_path.read_text()) if run_manifest_path.exists() else {}
        )
        checkpoint = (
            json.loads(checkpoint_path.read_text()) if checkpoint_path.exists() else {}
        )
        rows.append({
            "run": run.name,
            "status": run_manifest.get("status", "queued"),
            "step": int(checkpoint.get("step_number", 0)),
            "time": float(checkpoint.get("time", 0.0)),
            "git_sha": run_manifest.get("git_sha"),
            "checkpoint_age_seconds": (
                time.time() - checkpoint_path.stat().st_mtime
                if checkpoint_path.exists() else None
            ),
            "traceback": (run / "traceback.log").exists(),
        })
    started = [row for row in rows if row["status"] != "queued"]
    steps = [row["step"] for row in started]
    return {
        "campaign": str(campaign),
        "campaign_status": manifest.get("status", "missing"),
        "runs_total": len(rows),
        "runs_started": len(started),
        "status_counts": dict(sorted(Counter(row["status"] for row in rows).items())),
        "simulation_git_shas": sorted({
            row["git_sha"] for row in started if row["git_sha"]
        }),
        "checkpoint_step": {
            "minimum": min(steps) if steps else 0,
            "mean": sum(steps) / len(steps) if steps else 0.0,
            "maximum": max(steps) if steps else 0,
        },
        "oldest_checkpoint_age_seconds": max(
            (row["checkpoint_age_seconds"] for row in started
             if row["checkpoint_age_seconds"] is not None), default=None,
        ),
        "traceback_count": sum(row["traceback"] for row in rows),
        "runs": rows,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Report non-mutating campaign progress.")
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--compact", action="store_true")
    arguments = parser.parse_args()
    status = campaign_status(arguments.campaign)
    if arguments.compact:
        status.pop("runs")
    print(json.dumps(status, indent=2))
