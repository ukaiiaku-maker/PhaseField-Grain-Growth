from __future__ import annotations

import csv
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


class EventLedger:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._handle = self.path.open("a", newline="", encoding="utf-8")
        self._writer = csv.DictWriter(self._handle, fieldnames=EVENT_FIELDS, extrasaction="raise")
        if self.path.stat().st_size == 0:
            self._writer.writeheader()

    def write(self, record: dict[str, Any]) -> None:
        row = {name: record.get(name, "") for name in EVENT_FIELDS}
        self._writer.writerow(row)
        self._handle.flush()

    def close(self) -> None:
        self._handle.close()

    def __enter__(self) -> "EventLedger":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

