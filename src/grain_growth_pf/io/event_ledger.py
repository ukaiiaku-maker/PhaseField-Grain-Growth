from __future__ import annotations

import csv
import gzip
import json
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

_STRING_FIELDS = {
    "run_id", "event_id", "event_type", "grain_ids", "entity_id", "position",
    "barrier_type", "activation_volume", "burgers_vector_b", "Git_SHA",
}
_INTEGER_FIELDS = {"step", "seed", "hit_count", "required_hits_K"}


def event_ledger_path(run_dir: str | Path) -> Path:
    """Return the present Parquet, gzip, or plain event ledger for a run directory."""
    root = Path(run_dir)
    for name in ("events.parquet", "events.csv.gz", "events.csv"):
        candidate = root / name
        if candidate.exists():
            return candidate
    return root / "events.csv"


def event_ledger_has_rows(path: str | Path) -> bool:
    """Return whether an event ledger contains at least one committed data part."""
    path = Path(path)
    if not path.exists():
        return False
    if path.is_dir():
        return any(path.glob("part-*.parquet"))
    return path.stat().st_size > 0


def read_event_ledger(path: str | Path):
    """Load any supported event-ledger representation as a pandas DataFrame."""
    import pandas as pd

    path = Path(path)
    if path.is_dir() or path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path)


class EventLedger:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._handle = None
        self._writer = None
        self._columns: dict[str, list[Any]] | None = None
        self._part_count = 0
        self._open()

    @property
    def is_parquet(self) -> bool:
        return self.path.suffix == ".parquet"

    def _open(self) -> None:
        if self.is_parquet:
            self.path.mkdir(parents=True, exist_ok=True)
            parts = sorted(self.path.glob("part-*.parquet"))
            expected = [self.path / f"part-{index:08d}.parquet" for index in range(len(parts))]
            if parts != expected:
                raise ValueError(f"non-contiguous Parquet event-ledger parts in {self.path}")
            self._part_count = len(parts)
            self._columns = {name: [] for name in EVENT_FIELDS}
            return
        empty = not self.path.exists() or self.path.stat().st_size == 0
        if self.path.suffix == ".gz":
            self._handle = gzip.open(
                self.path, "at", newline="", encoding="utf-8", compresslevel=1
            )
        else:
            self._handle = self.path.open(
                "a", newline="", encoding="utf-8", buffering=1024 * 1024
            )
        self._writer = csv.DictWriter(self._handle, fieldnames=EVENT_FIELDS, extrasaction="raise")
        if empty:
            self._writer.writeheader()

    def write(self, record: dict[str, Any]) -> None:
        if self.is_parquet:
            assert self._columns is not None
            for name in EVENT_FIELDS:
                value = record.get(name)
                if isinstance(value, str) and value == "":
                    value = None
                elif name in _STRING_FIELDS and value is not None and not isinstance(value, str):
                    serializable = value.tolist() if hasattr(value, "tolist") else value
                    try:
                        value = json.dumps(serializable)
                    except TypeError:
                        value = str(value)
                self._columns[name].append(value)
            return
        row = {name: record.get(name, "") for name in EVENT_FIELDS}
        self._writer.writerow(row)

    @staticmethod
    def _parquet_schema():
        import pyarrow as pa

        return pa.schema([
            pa.field(
                name,
                pa.string() if name in _STRING_FIELDS
                else pa.int64() if name in _INTEGER_FIELDS
                else pa.float64(),
            )
            for name in EVENT_FIELDS
        ])

    def _checkpoint_parquet(self) -> int:
        assert self._columns is not None
        rows = len(self._columns[EVENT_FIELDS[0]])
        if not rows:
            return self._part_count
        import pyarrow as pa
        import pyarrow.parquet as pq

        table = pa.Table.from_pydict(self._columns, schema=self._parquet_schema())
        final_path = self.path / f"part-{self._part_count:08d}.parquet"
        temporary_path = self.path / f".{final_path.name}.{os.getpid()}.tmp"
        pq.write_table(table, temporary_path, compression="zstd", compression_level=1)
        with temporary_path.open("rb") as handle:
            os.fsync(handle.fileno())
        os.replace(temporary_path, final_path)
        directory_fd = os.open(self.path, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
        self._part_count += 1
        self._columns = {name: [] for name in EVENT_FIELDS}
        return self._part_count

    def checkpoint(self) -> int:
        """Commit one restart-consistent generation and return its storage extent."""
        if self.is_parquet:
            return self._checkpoint_parquet()
        self.close()
        with self.path.open("rb") as handle:
            os.fsync(handle.fileno())
        offset = self.path.stat().st_size
        self._open()
        return offset

    def truncate(self, offset: int) -> None:
        """Discard rows or compressed members newer than an authoritative checkpoint."""
        if self.is_parquet:
            parts = sorted(self.path.glob("part-*.parquet"))
            if offset < 0 or offset > len(parts):
                raise ValueError(
                    f"invalid event-ledger checkpoint part {offset} for {len(parts)} parts"
                )
            for part in parts[offset:]:
                part.unlink()
            for orphan in self.path.glob(".part-*.tmp"):
                orphan.unlink()
            self._part_count = offset
            self._columns = {name: [] for name in EVENT_FIELDS}
            directory_fd = os.open(self.path, os.O_RDONLY)
            try:
                os.fsync(directory_fd)
            finally:
                os.close(directory_fd)
            return
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
        if self.is_parquet:
            self._checkpoint_parquet()
            return
        if self._handle is not None and not self._handle.closed:
            self._handle.close()

    def __enter__(self) -> "EventLedger":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()
