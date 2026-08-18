#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from grain_growth_pf.disconnections.mode import K_B_EV
from grain_growth_pf.io.event_ledger import event_ledger_path, iter_event_ledger


@dataclass
class Stats:
    count: int = 0
    barrierless: int = 0
    total: float = 0.0
    minimum: float = np.inf
    maximum: float = -np.inf
    sample_values: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    sample_keys: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))

    def add(self, values: np.ndarray, tol: float, rng: np.random.Generator, sample_size: int) -> None:
        values = np.asarray(values, dtype=float)
        if values.size == 0:
            return
        self.count += int(values.size)
        self.barrierless += int(np.count_nonzero(values <= tol))
        self.total += float(values.sum())
        self.minimum = min(self.minimum, float(values.min()))
        self.maximum = max(self.maximum, float(values.max()))
        keys = rng.random(values.size)
        if self.sample_values.size:
            values = np.concatenate((self.sample_values, values))
            keys = np.concatenate((self.sample_keys, keys))
        if values.size > sample_size:
            keep = np.argpartition(keys, sample_size - 1)[:sample_size]
            values = values[keep]
            keys = keys[keep]
        self.sample_values = values
        self.sample_keys = keys


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stream event ledgers and quantify the attempt-limited (effective barrier = 0) fraction."
    )
    parser.add_argument("campaign")
    parser.add_argument("--regime", action="append", dest="regimes")
    parser.add_argument("--event-type", action="append", dest="event_types")
    parser.add_argument("--output", required=True)
    parser.add_argument("--tolerance-ev", type=float, default=1.0e-12)
    parser.add_argument("--sample-size", type=int, default=100_000)
    args = parser.parse_args()

    campaign = Path(args.campaign)
    manifest = json.loads((campaign / "campaign_manifest.json").read_text())
    selected_regimes = set(args.regimes or [])
    selected_event_types = set(args.event_types or [])
    accum: dict[tuple[str, float], Stats] = {}
    rngs: dict[tuple[str, float], np.random.Generator] = {}
    runs_used: dict[tuple[str, float], int] = {}

    for raw_run in manifest.get("runs", []):
        run = Path(raw_run)
        run_manifest_path = run / "manifest.json"
        if not run_manifest_path.exists():
            continue
        run_manifest = json.loads(run_manifest_path.read_text())
        if run_manifest.get("status") != "completed":
            continue
        config = run_manifest["config"]
        regime = str(config["regime"])
        temperature = float(config["pf"]["temperature"])
        if selected_regimes and regime not in selected_regimes:
            continue
        key = (regime, temperature)
        if key not in accum:
            accum[key] = Stats()
            seed = (sum(ord(ch) for ch in regime) * 1009 + int(round(temperature * 10))) % (2**32)
            rngs[key] = np.random.default_rng(seed)
            runs_used[key] = 0
        runs_used[key] += 1
        ledger = event_ledger_path(run)
        if not ledger.exists():
            continue
        for frame in iter_event_ledger(
            ledger, columns=["event_type", "DeltaG0", "effective_DeltaG"], batch_size=250_000
        ):
            if selected_event_types:
                frame = frame[frame["event_type"].isin(selected_event_types)]
                if frame.empty:
                    continue
            bare = pd.to_numeric(frame["DeltaG0"], errors="coerce")
            effective = pd.to_numeric(frame["effective_DeltaG"], errors="coerce")
            valid = bare.notna() & effective.notna()
            values = effective.loc[valid].to_numpy(dtype=float)
            accum[key].add(values, args.tolerance_ev, rngs[key], args.sample_size)

    rows = []
    for (regime, temperature), stats in sorted(accum.items()):
        sample = stats.sample_values
        quantiles = np.quantile(sample, [0.1, 0.5, 0.9]) if sample.size else [np.nan] * 3
        kbt = K_B_EV * temperature
        rows.append({
            "regime": regime,
            "temperature_K": temperature,
            "completed_runs": runs_used[(regime, temperature)],
            "barrier_events": stats.count,
            "barrierless_events": stats.barrierless,
            "barrierless_fraction": stats.barrierless / stats.count if stats.count else np.nan,
            "effective_DeltaG_mean_eV": stats.total / stats.count if stats.count else np.nan,
            "effective_DeltaG_min_eV": stats.minimum if stats.count else np.nan,
            "effective_DeltaG_q10_eV": quantiles[0],
            "effective_DeltaG_q50_eV": quantiles[1],
            "effective_DeltaG_q90_eV": quantiles[2],
            "effective_DeltaG_max_eV": stats.maximum if stats.count else np.nan,
            "median_barrier_over_kBT": quantiles[1] / kbt if sample.size else np.nan,
            "sampled_for_quantiles": int(sample.size),
        })

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result = pd.DataFrame(rows)
    result.to_csv(output, index=False)
    if len(result):
        print(result.to_string(index=False))
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
