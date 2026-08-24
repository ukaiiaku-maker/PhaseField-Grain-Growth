from __future__ import annotations

import csv
from collections import deque
import json
import os
from pathlib import Path
from typing import Any, Iterable, Mapping


TRACE_FIELDS = (
    "run_id", "step", "time", "temperature", "seed", "event_ids",
    "entity_type", "entity_id", "grain_ids", "position",
    "gb_arclength_coordinate", "adjoining_tj_ids", "local_normal_velocity",
    "local_signed_normal_displacement", "gb_length", "domain_length",
    "curvature", "G_pending", "T_pending", "C_pending", "climb_stage",
    "shear_state_s", "tau_int", "free_volume_signed_inventory",
    "p_cap", "p_chem", "p_shear", "p_event", "p_net", "p_applied_total", "chi_s",
    "local_shear_energy",
    "p_cap_V_n", "tau_V_tau", "Delta_mu_v_N_v", "W_total", "DeltaG0",
    "DeltaG_eff", "instantaneous_rate", "cumulative_hazard",
    "hazard_threshold", "N_required", "N_accommodated_GB",
    "N_accommodated_TJ", "N_stored", "conservation_residual", "tj_velocity",
    "tj_compatibility_residual", "residual_burgers", "tj_sink_admissible",
    "signed_point_defect_flux",
)

EVENT_INDEX_FIELDS = (
    "run_id", "event_id", "event_type", "sink_path", "event_step",
    "event_time", "entity_type", "entity_id", "grain_ids", "position",
    "gb_arclength_coordinate", "adjoining_tj_ids", "trace_start_step",
    "trace_end_step", "trace_stride", "neighborhood_entities",
    "instantaneous_rate", "cumulative_hazard", "hazard_threshold", "DeltaG0",
    "DeltaG_eff", "signed_point_defect_flux", "tj_velocity",
    "tj_compatibility_residual", "residual_burgers", "tj_sink_admissible",
)

TRACE_STRING_FIELDS = {
    "run_id", "event_ids", "entity_type", "entity_id", "grain_ids", "position",
    "adjoining_tj_ids", "climb_stage", "tj_velocity",
    "tj_compatibility_residual", "residual_burgers", "tj_sink_admissible",
}
EVENT_STRING_FIELDS = {
    "run_id", "event_id", "event_type", "sink_path", "entity_type", "entity_id",
    "grain_ids", "position", "adjoining_tj_ids", "neighborhood_entities",
    "tj_velocity", "tj_compatibility_residual", "residual_burgers",
    "tj_sink_admissible",
}


def trace_entity_key(entity_type: str, entity_id: str) -> str:
    return f"{entity_type}::{entity_id}"


class EventTraceRecorder:
    """Merge event windows into unique solver-step/local-entity trace rows.

    The event index stores every trigger and its neighborhood.  The state table
    stores each ``(entity_type, entity_id, step)`` at most once, so overlapping
    pre/post windows do not duplicate local state.
    """

    def __init__(
        self,
        output_dir: str | Path,
        *,
        run_id: str,
        pre_steps: int,
        post_steps: int,
        stride: int,
        output_format: str = "csv",
        resume: bool = False,
        checkpoint_state: Mapping[str, Any] | None = None,
    ) -> None:
        if pre_steps < 0 or post_steps < 0 or stride <= 0:
            raise ValueError("event trace pre/post steps must be nonnegative and stride positive")
        self.run_id = str(run_id)
        self.pre_steps = int(pre_steps)
        self.post_steps = int(post_steps)
        self.stride = int(stride)
        self.output_format = str(output_format).lower()
        if self.output_format not in {"csv", "parquet"}:
            raise ValueError("event trace format must be csv or parquet")
        self.output_dir = Path(output_dir)
        suffix = "parquet" if self.output_format == "parquet" else "csv"
        self.state_path = self.output_dir / f"event_traces.{suffix}"
        self.event_path = self.output_dir / f"event_trace_events.{suffix}"
        state = dict(checkpoint_state or {})
        if resume and not state:
            raise ValueError("instrumented restart lacks event-trace checkpoint state")
        if resume:
            if str(state.get("output_format", "csv")) != self.output_format:
                raise ValueError("event-trace format must match the checkpoint")
            restored_shape = (
                int(state.get("pre_steps", -1)),
                int(state.get("post_steps", -1)),
                int(state.get("stride", -1)),
            )
            requested_shape = (self.pre_steps, self.post_steps, self.stride)
            if restored_shape != requested_shape:
                raise ValueError(
                    "event-trace window settings must match the checkpoint: "
                    f"checkpoint={restored_shape}, requested={requested_shape}"
                )
        self._state_handle = self._event_handle = None
        self._state_writer = self._event_writer = None
        self._state_rows: list[dict[str, Any]] = []
        self._event_rows: list[dict[str, Any]] = []
        self._state_part_count = int(state.get("state_part_count", 0))
        self._event_part_count = int(state.get("event_part_count", 0))
        if self.output_format == "parquet":
            self._restore_parquet_parts(
                self.state_path, resume, self._state_part_count
            )
            self._restore_parquet_parts(
                self.event_path, resume, self._event_part_count
            )
        else:
            self._state_handle, self._state_writer = self._open_stream(
                self.state_path, TRACE_FIELDS, resume, state.get("state_offset")
            )
            self._event_handle, self._event_writer = self._open_stream(
                self.event_path, EVENT_INDEX_FIELDS, resume, state.get("event_offset")
            )
        self._buffer: deque[tuple[int, dict[str, dict[str, Any]]]] = deque(
            (int(item["step"]), dict(item["rows"]))
            for item in state.get("buffer", [])
        )
        self._active: dict[str, dict[str, int]] = {
            str(key): {str(event_id): int(end) for event_id, end in value.items()}
            for key, value in state.get("active", {}).items()
        }
        self._pending: list[dict[str, Any]] = [
            dict(item) for item in state.get("pending", [])
        ]
        self._recent_written: set[tuple[str, int]] = {
            (str(key), int(step)) for key, step in state.get("recent_written", [])
        }

    @staticmethod
    def _restore_parquet_parts(path: Path, resume: bool, part_count: int) -> None:
        if not resume:
            path.mkdir(parents=True, exist_ok=False)
            return
        if not path.is_dir():
            raise ValueError(f"missing restart-safe event trace dataset {path}")
        parts = sorted(path.glob("part-*.parquet"))
        expected = [path / f"part-{index:08d}.parquet" for index in range(len(parts))]
        if parts != expected or part_count < 0 or part_count > len(parts):
            raise ValueError(f"invalid event trace Parquet generation under {path}")
        for orphan in parts[part_count:]:
            orphan.unlink()

    @staticmethod
    def _open_stream(
        path: Path,
        fields: tuple[str, ...],
        resume: bool,
        offset: Any,
    ) -> tuple[Any, csv.DictWriter]:
        if resume:
            if offset is None or not path.exists():
                raise ValueError(f"missing restart-safe event trace stream {path}")
            with path.open("r+b") as raw:
                size = raw.seek(0, os.SEEK_END)
                position = int(offset)
                if position < 0 or position > size:
                    raise ValueError(f"invalid event trace offset {position} for {size} bytes")
                raw.truncate(position)
                raw.flush()
                os.fsync(raw.fileno())
            handle = path.open("a", newline="", encoding="utf-8")
        else:
            handle = path.open("w", newline="", encoding="utf-8")
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        if not resume:
            writer.writeheader()
            handle.flush()
        return handle, writer

    @property
    def requested_keys(self) -> set[str]:
        keys = set(self._active)
        for event in self._pending:
            keys.update(map(str, event["neighborhood_keys"]))
        return keys

    def trigger(self, record: Mapping[str, Any], neighborhood_keys: Iterable[str]) -> None:
        event_id = str(record.get("event_id", ""))
        if not event_id:
            raise ValueError("event trace trigger requires an event_id")
        step = int(record["step"])
        keys = sorted(set(map(str, neighborhood_keys)))
        self._pending.append({
            "record": dict(record),
            "event_id": event_id,
            "step": step,
            "start_step": step - self.pre_steps,
            "end_step": step + self.post_steps,
            "neighborhood_keys": keys,
        })

    def capture(self, step: int, rows_by_key: Mapping[str, Mapping[str, Any]]) -> None:
        current = int(step)
        sample_due = current % self.stride == 0 or bool(self._pending)
        rows = {str(key): dict(value) for key, value in rows_by_key.items()}
        if sample_due:
            self._buffer.append((current, rows))

        for event in self._pending:
            record = dict(event["record"])
            key_list = list(map(str, event["neighborhood_keys"]))
            self._write_event_index({
                "run_id": self.run_id,
                "event_id": event["event_id"],
                "event_type": record.get("event_type", ""),
                "sink_path": record.get("sink_path", ""),
                "event_step": event["step"],
                "event_time": record.get("time", ""),
                "entity_type": "TJ" if str(record.get("entity_id", "")).startswith("tj:") else "GB",
                "entity_id": record.get("entity_id", ""),
                "grain_ids": record.get("grain_ids", ""),
                "position": self._json(record.get("position", "")),
                "gb_arclength_coordinate": record.get("gb_arclength_coordinate", ""),
                "adjoining_tj_ids": self._json(record.get("adjoining_tj_ids", "")),
                "trace_start_step": event["start_step"],
                "trace_end_step": event["end_step"],
                "trace_stride": self.stride,
                "neighborhood_entities": json.dumps(key_list, separators=(",", ":")),
                "instantaneous_rate": record.get("instantaneous_rate", ""),
                "cumulative_hazard": record.get("cumulative_hazard", ""),
                "hazard_threshold": record.get("random_hazard_threshold", ""),
                "DeltaG0": record.get("DeltaG0", ""),
                "DeltaG_eff": record.get("effective_DeltaG", ""),
                "signed_point_defect_flux": record.get("signed_defect_quota", record.get("Nv", "")),
                "tj_velocity": self._json(record.get("tj_velocity_compatible", "")),
                "tj_compatibility_residual": self._json(record.get("compatibility_residual", "")),
                "residual_burgers": self._json(record.get("burgers_after", record.get("burgers_vector_b", ""))),
                "tj_sink_admissible": record.get("candidate_allowed", ""),
            })
            for buffered_step, buffered_rows in self._buffer:
                if event["start_step"] <= buffered_step <= current:
                    for key in key_list:
                        row = buffered_rows.get(key)
                        if row is not None:
                            self._write_state(key, buffered_step, row, [event["event_id"]])
            for key in key_list:
                self._active.setdefault(key, {})[event["event_id"]] = int(event["end_step"])
        self._pending.clear()

        if sample_due:
            for key, events in list(self._active.items()):
                row = rows.get(key)
                if row is not None:
                    event_ids = sorted(event_id for event_id, end in events.items() if current <= end)
                    if event_ids:
                        self._write_state(key, current, row, event_ids)
                retained = {event_id: end for event_id, end in events.items() if current < end}
                if retained:
                    self._active[key] = retained
                else:
                    self._active.pop(key, None)

        minimum = current - self.pre_steps
        while self._buffer and self._buffer[0][0] < minimum:
            self._buffer.popleft()
        self._recent_written = {
            (key, recorded_step)
            for key, recorded_step in self._recent_written
            if recorded_step >= minimum
        }

    def _write_event_index(self, row: Mapping[str, Any]) -> None:
        output = {name: row.get(name, "") for name in EVENT_INDEX_FIELDS}
        if self.output_format == "parquet":
            self._event_rows.append(output)
            if len(self._event_rows) >= 25_000:
                self._flush_parquet("event")
        else:
            assert self._event_writer is not None
            self._event_writer.writerow(output)

    def _write_state(
        self,
        key: str,
        step: int,
        row: Mapping[str, Any],
        event_ids: Iterable[str],
    ) -> None:
        identity = (key, int(step))
        if identity in self._recent_written:
            return
        output = {name: row.get(name, "") for name in TRACE_FIELDS}
        output["run_id"] = self.run_id
        output["step"] = int(step)
        output["event_ids"] = ";".join(sorted(set(map(str, event_ids))))
        if self.output_format == "parquet":
            self._state_rows.append(output)
            if len(self._state_rows) >= 25_000:
                self._flush_parquet("state")
        else:
            assert self._state_writer is not None
            self._state_writer.writerow(output)
        self._recent_written.add(identity)

    def _flush_parquet(self, kind: str) -> None:
        import pyarrow as pa
        import pyarrow.parquet as pq

        is_state = kind == "state"
        rows = self._state_rows if is_state else self._event_rows
        if not rows:
            return
        fields = TRACE_FIELDS if is_state else EVENT_INDEX_FIELDS
        strings = TRACE_STRING_FIELDS if is_state else EVENT_STRING_FIELDS
        integer_fields = {"step", "seed"} if is_state else {
            "event_step", "trace_start_step", "trace_end_step", "trace_stride"
        }
        schema = pa.schema([
            pa.field(name, pa.string() if name in strings else (
                pa.int64() if name in integer_fields else pa.float64()
            ))
            for name in fields
        ])
        normalized = []
        for row in rows:
            item: dict[str, Any] = {}
            for name in fields:
                value = row.get(name)
                blank = value is None or (isinstance(value, str) and value == "")
                if blank:
                    item[name] = None
                elif name in strings:
                    item[name] = str(value)
                elif name in integer_fields:
                    item[name] = int(value)
                else:
                    item[name] = float(value)
            normalized.append(item)
        part_count = self._state_part_count if is_state else self._event_part_count
        path = (self.state_path if is_state else self.event_path) / f"part-{part_count:08d}.parquet"
        pq.write_table(
            pa.Table.from_pylist(normalized, schema=schema), path,
            compression="zstd", use_dictionary=True,
        )
        rows.clear()
        if is_state:
            self._state_part_count += 1
        else:
            self._event_part_count += 1

    @staticmethod
    def _json(value: Any) -> str:
        if value is None or (isinstance(value, str) and value == ""):
            return ""
        if isinstance(value, str):
            return value
        serializable = value.tolist() if hasattr(value, "tolist") else value
        return json.dumps(serializable, separators=(",", ":"))

    def checkpoint_state(self) -> dict[str, Any]:
        if self.output_format == "parquet":
            self._flush_parquet("state")
            self._flush_parquet("event")
        else:
            for handle in (self._state_handle, self._event_handle):
                assert handle is not None
                handle.flush()
                os.fsync(handle.fileno())
        state = {
            "pre_steps": self.pre_steps,
            "post_steps": self.post_steps,
            "stride": self.stride,
            "output_format": self.output_format,
            "buffer": [
                {"step": step, "rows": rows} for step, rows in self._buffer
            ],
            "active": self._active,
            "pending": self._pending,
            "recent_written": sorted([key, step] for key, step in self._recent_written),
        }
        if self.output_format == "parquet":
            state.update({
                "state_part_count": self._state_part_count,
                "event_part_count": self._event_part_count,
            })
        else:
            assert self._state_handle is not None and self._event_handle is not None
            state.update({
                "state_offset": self._state_handle.tell(),
                "event_offset": self._event_handle.tell(),
            })
        return state

    def close(self) -> None:
        if self.output_format == "parquet":
            self._flush_parquet("state")
            self._flush_parquet("event")
            return
        for handle in (self._state_handle, self._event_handle):
            if handle is not None and not handle.closed:
                handle.flush()
                handle.close()
