from __future__ import annotations

import csv
import gzip
import os
from pathlib import Path
from typing import Any, Iterable


EVENT_FIELDS = (
    "run_id", "time", "step", "temperature", "seed", "event_id", "event_type",
    "grain_ids", "entity_id", "position", "geometry_measure_Q", "grain_size",
    "curvature", "local_velocity", "barrier_type", "DeltaG0", "effective_DeltaG",
    "activation_volume", "local_shear_stress", "local_normal_free_volume_stress",
    "shear_state_s", "free_volume_state_q", "Ns", "nu0", "instantaneous_rate",
    "cumulative_hazard", "random_hazard_threshold", "hit_count", "required_hits_K",
    "release_Delta_s", "release_Delta_q", "GB_area_change", "TJ_travel",
    "point_defect_quota", "normal_step_h", "burgers_vector_b", "Nv",
    "shear_strain_increment", "volumetric_strain_increment", "packet_size", "Git_SHA",
)


def event_ledger_path(run_dir: str | Path) -> Path:
    """Return the present plain or gzip event ledger for a run directory."""
    root = Path(run_dir)
    for name in ("events.csv.gz", "events.csv"):
        candidate = root / name
        if candidate.exists():
            return candidate
    return root / "events.csv"


class EventLedger:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._handle = None
        self._writer = None
        self._open()

    def _open(self) -> None:
        empty = not self.path.exists() or self.path.stat().st_size == 0
        if self.path.suffix == ".gz":
            self._handle = gzip.open(
                self.path, "at", newline="", encoding="utf-8", compresslevel=6
            )
        else:
            self._handle = self.path.open(
                "a", newline="", encoding="utf-8", buffering=1024 * 1024
            )
        self._writer = csv.DictWriter(self._handle, fieldnames=EVENT_FIELDS, extrasaction="raise")
        if empty:
            self._writer.writeheader()

    def write(self, record: dict[str, Any]) -> None:
        row = {name: record.get(name, "") for name in EVENT_FIELDS}
        self._writer.writerow(row)

    def checkpoint(self) -> int:
        """Commit one restart-consistent ledger generation and return its byte extent."""
        self.close()
        with self.path.open("rb") as handle:
            os.fsync(handle.fileno())
        offset = self.path.stat().st_size
        self._open()
        return offset

    def truncate(self, offset: int) -> None:
        """Discard rows or compressed members newer than an authoritative checkpoint."""
        self.close()
        size = self.path.stat().st_size
        if offset < 0 or offset > size:
            raise ValueError(f"invalid event-ledger checkpoint offset {offset} for {size} bytes")
        with self.path.open("r+b") as handle:
            handle.truncate(offset)
            handle.flush()
            os.fsync(handle.fileno())
        self._open()

    def close(self) -> None:
        if self._handle is not None and not self._handle.closed:
            self._handle.close()

    def __enter__(self) -> "EventLedger":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()
