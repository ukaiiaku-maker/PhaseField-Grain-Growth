#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

import pandas as pd


def _root_manifest(root: Path) -> dict[str, Any]:
    for name in ("campaign_manifest.json", "video_manifest.json"):
        path = root / name
        if path.exists():
            return json.loads(path.read_text())
    return {}


def _run_paths(root: Path, root_manifest: dict[str, Any]) -> list[Path]:
    if root_manifest.get("runs"):
        return [
            Path(item["path"] if isinstance(item, dict) else item)
            for item in root_manifest["runs"]
        ]
    return sorted(path.parent for path in root.glob("*/manifest.json"))


def _movie(run: Path) -> Path | None:
    candidates = sorted(run.glob("final.*")) + sorted(run.glob("diagnostic.*"))
    return next((path for path in candidates if path.suffix.lower() in {".mp4", ".gif"}), None)


def _restart(run: Path, manifest: dict[str, Any], root_manifest: dict[str, Any]) -> dict[str, Any]:
    if manifest.get("restart_provenance"):
        return dict(manifest["restart_provenance"])
    provenance = root_manifest.get("restart_provenance", {})
    return dict(provenance.get(str(run), {}))


def _status(manifest: dict[str, Any], final_grains: int, ending_step: int) -> str:
    config = manifest["config"]
    target = int(config.get("termination_grains", 0))
    ceiling = int(config.get("max_steps", 0))
    if str(manifest.get("status", "")) != "completed":
        return str(manifest.get("status", "incomplete"))
    if target and final_grains <= target:
        return "target_reached"
    if ceiling and ending_step >= ceiling:
        return "censored_step_ceiling"
    return "completed_before_ceiling"


def summarize(roots: list[Path]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for root in roots:
        root_manifest = _root_manifest(root)
        for run in _run_paths(root, root_manifest):
            manifest_path = run / "manifest.json"
            checkpoint_path = run / "checkpoint.json"
            if not manifest_path.exists() or not checkpoint_path.exists():
                continue
            manifest = json.loads(manifest_path.read_text())
            checkpoint = json.loads(checkpoint_path.read_text())
            config = manifest["config"]
            restart = _restart(run, manifest, root_manifest)
            ending_step = int(checkpoint["step_number"])
            final_grains = int(manifest.get("final_grains", -1))
            movie = _movie(run)
            rows.append({
                "campaign": root.name,
                "regime": str(config["regime"]),
                "path": str(run),
                "source_sha": str(manifest.get("git_sha", root_manifest.get("git_sha", ""))),
                "starting_step": int(restart.get("checkpoint_step", 0)),
                "starting_time": float(restart.get("checkpoint_time", 0.0)),
                "starting_grains": int(restart.get(
                    "starting_grains", manifest.get("grains_after_equilibration", 200)
                )),
                "ending_step": ending_step,
                "ending_time": float(checkpoint["time"]),
                "ending_grains": final_grains,
                "completion_status": _status(manifest, final_grains, ending_step),
                "frame_count": len(list((run / "frames").glob("frame-*.npz"))),
                "movie_path": str(movie) if movie else "MISSING",
            })
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(["campaign", "regime"]).reset_index(drop=True)


def markdown(table: pd.DataFrame, document: Path) -> str:
    lines = ["# Area-loss climb execution and movie index", ""]
    for campaign, frame in table.groupby("campaign", sort=True):
        lines.extend([f"## {campaign}", "", "| Regime | Start | End | Grains | Status | Frames | Movie |", "|---|---:|---:|---:|---|---:|---|"])
        for row in frame.itertuples():
            if row.movie_path == "MISSING":
                movie = "MISSING"
            else:
                relative = os.path.relpath(Path(row.movie_path), start=document.parent)
                movie = f"[{Path(row.movie_path).name}]({relative})"
            lines.append(
                f"| {row.regime} | {row.starting_step} / {row.starting_time:.6g} | "
                f"{row.ending_step} / {row.ending_time:.6g} | "
                f"{row.starting_grains}→{row.ending_grains} | {row.completion_status} | "
                f"{row.frame_count} | {movie} |"
            )
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="Index area-loss execution provenance and movies.")
    parser.add_argument("campaign", nargs="+", type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--markdown", required=True, type=Path)
    args = parser.parse_args()
    table = summarize(args.campaign)
    if table.empty:
        raise SystemExit("no run manifests found")
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.csv, index=False)
    args.markdown.write_text(markdown(table, args.markdown) + "\n")
    print(args.csv)
    print(args.markdown)


if __name__ == "__main__":
    main()
