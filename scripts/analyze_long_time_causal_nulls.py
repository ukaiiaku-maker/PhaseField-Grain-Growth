#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from grain_growth_pf.analysis.causal import trace_center_causal_null


WINDOWS = (
    (1.00, 1.20), (1.15, 1.40), (1.30, 1.60), (1.50, 1.85),
    (1.75, 2.15), (2.00, 2.50), (2.35, 2.90), (2.70, 3.20),
)


def _family(event_type: str) -> str:
    if event_type == "compatibility_release":
        return "G_compatibility_release"
    if event_type == "tj_compatibility_release":
        return "T_compatibility_release"
    if event_type == "gb_sink_completion":
        return "C_GB_sink_completion"
    if event_type == "tj_sink_completion":
        return "C_TJ_sink_completion"
    return event_type


def _stable_seed(*parts: object) -> int:
    digest = hashlib.sha256("|".join(map(str, parts)).encode()).digest()
    return int.from_bytes(digest[:8], "little")


def _run_rows(run: Path, trajectory: pd.DataFrame, shuffles: int) -> list[dict]:
    event_root = run / "event_trace_events.parquet"
    trace_root = run / "event_traces.parquet"
    if not list(event_root.glob("*.parquet")) or not list(trace_root.glob("*.parquet")):
        return []
    events = pd.read_parquet(
        event_root,
        columns=[
            "event_id", "event_step", "event_type", "trace_start_step",
            "trace_end_step",
        ],
    )
    records: dict[str, list[tuple[int, int, int, str]]] = defaultdict(list)
    for event in events.itertuples(index=False):
        records[str(event.event_id)].append(
            (
                int(event.event_step), int(event.trace_start_step),
                int(event.trace_end_step), _family(str(event.event_type)),
            )
        )
    known_steps = trajectory.step.to_numpy(int)
    known_x = trajectory.G_population_over_G0.to_numpy(float)
    profile_sums: dict[tuple[int, str, str], dict[int, list[float]]] = defaultdict(
        lambda: defaultdict(lambda: [0.0, 0.0])
    )
    for part in sorted(trace_root.glob("*.parquet")):
        parquet = pq.ParquetFile(part)
        for batch in parquet.iter_batches(
            batch_size=250_000,
            columns=["step", "event_ids", "entity_type", "local_normal_velocity"],
        ):
            frame = batch.to_pandas()
            for trace in frame.itertuples(index=False):
                if not isinstance(trace.event_ids, str) or not trace.event_ids:
                    continue
                velocity = float(trace.local_normal_velocity)
                if not np.isfinite(velocity):
                    continue
                for event_id in trace.event_ids.split(";"):
                    candidates = [
                        item for item in records.get(event_id, [])
                        if item[1] <= int(trace.step) <= item[2]
                    ]
                    if not candidates:
                        continue
                    event_step, _, _, family = min(
                        candidates,
                        key=lambda item: (abs(int(trace.step) - item[0]), item[0]),
                    )
                    key = (event_step, family, str(trace.entity_type))
                    relative = int(trace.step) - event_step
                    profile_sums[key][relative][0] += abs(velocity)
                    profile_sums[key][relative][1] += 1.0
    grouped: dict[tuple[float, float, str, str], list[dict[int, float]]] = defaultdict(list)
    for (event_step, family, entity_type), values in profile_sums.items():
        x = float(np.interp(event_step, known_steps, known_x))
        window = next(((lo, hi) for lo, hi in WINDOWS if lo <= x <= hi), None)
        if window is None:
            continue
        profile = {
            relative: total / count
            for relative, (total, count) in values.items()
            if count
        }
        grouped[(window[0], window[1], family, entity_type)].append(profile)
    rows = []
    regime = str(trajectory.regime.iloc[0])
    for (lo, hi, family, entity_type), profiles in sorted(grouped.items()):
        result = trace_center_causal_null(
            profiles,
            shuffles=shuffles,
            seed=_stable_seed(regime, lo, hi, family, entity_type),
        )
        rows.append(
            {
                "regime": regime,
                "window_start_x": lo,
                "window_end_x": hi,
                "event_family": family,
                "observed_entity_type": entity_type,
                **result,
                "null_definition": "within-event-trace pseudo-center outside +/-5 steps",
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--output")
    parser.add_argument("--shuffles", type=int, default=200)
    args = parser.parse_args()
    root = Path(args.campaign)
    output = Path(args.output) if args.output else root / "analysis" / "causal_nulls_vs_grain_size.csv"
    trajectory = pd.read_csv(root / "analysis" / "grain_size_trajectory.csv")
    manifest = json.loads((root / "campaign_manifest.json").read_text(encoding="utf-8"))
    rows = []
    for raw in manifest["runs"]:
        run = Path(raw)
        regime = json.loads((run / "manifest.json").read_text(encoding="utf-8"))["config"]["regime"]
        rows.extend(_run_rows(run, trajectory[trajectory.regime == regime], args.shuffles))
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output, index=False)
    print(output)


if __name__ == "__main__":
    main()
