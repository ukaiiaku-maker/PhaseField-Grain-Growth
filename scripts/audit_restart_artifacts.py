#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from grain_growth_pf.io.provenance import file_sha256, git_sha


def audit_campaign(campaign_dir: str | Path) -> int:
    campaign_dir = Path(campaign_dir)
    campaign_path = campaign_dir / "campaign_manifest.json"
    campaign = json.loads(campaign_path.read_text())
    audited = 0
    for raw_run in campaign["runs"]:
        run = Path(raw_run)
        manifest_path = run / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        artifacts = []
        for name in ("checkpoint.npz", "checkpoint.json"):
            path = run / name
            if path.exists():
                artifacts.append({
                    "path": str(path),
                    "sha256": file_sha256(path),
                    "size_bytes": path.stat().st_size,
                })
        manifest["restart_artifacts"] = artifacts
        manifest["artifact_audit_git_sha"] = git_sha()
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        audited += 1
    campaign["restart_artifact_audit"] = {
        "git_sha": git_sha(), "runs_audited": audited,
    }
    campaign_path.write_text(json.dumps(campaign, indent=2) + "\n")
    return audited


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add SHA-256 restart-file provenance to an existing immutable campaign."
    )
    parser.add_argument("campaign")
    arguments = parser.parse_args()
    print(audit_campaign(arguments.campaign))
