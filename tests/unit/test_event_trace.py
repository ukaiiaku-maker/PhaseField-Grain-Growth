from __future__ import annotations

import csv

import pytest

from grain_growth_pf.io.event_trace import EventTraceRecorder, trace_entity_key


def _row(entity_id: str, value: float) -> dict[str, object]:
    return {
        "time": value,
        "entity_type": "GB",
        "entity_id": entity_id,
        "grain_ids": "1;2",
        "local_normal_velocity": value,
    }


def _read(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_event_trace_merges_overlapping_windows_without_duplicate_state_rows(tmp_path):
    key = trace_entity_key("GB", "gb:1-2:0")
    recorder = EventTraceRecorder(
        tmp_path, run_id="run", pre_steps=2, post_steps=2, stride=1
    )
    for step in (1, 2):
        recorder.capture(step, {key: _row("gb:1-2:0", float(step))})
    recorder.trigger(
        {
            "event_id": "event-1",
            "event_type": "compatibility_release",
            "entity_id": "gb:1-2:0",
            "step": 3,
            "time": 3.0,
        },
        [key],
    )
    recorder.capture(3, {key: _row("gb:1-2:0", 3.0)})
    recorder.trigger(
        {
            "event_id": "event-2",
            "event_type": "gb_sink_completion",
            "entity_id": "gb:1-2:0",
            "step": 4,
            "time": 4.0,
        },
        [key],
    )
    for step in (4, 5, 6):
        recorder.capture(step, {key: _row("gb:1-2:0", float(step))})
    recorder.close()

    states = _read(tmp_path / "event_traces.csv")
    events = _read(tmp_path / "event_trace_events.csv")
    identities = [(row["entity_id"], int(row["step"])) for row in states]
    assert identities == [("gb:1-2:0", step) for step in range(1, 7)]
    assert len(identities) == len(set(identities))
    assert [row["event_id"] for row in events] == ["event-1", "event-2"]
    assert events[0]["trace_start_step"] == "1"
    assert events[1]["trace_end_step"] == "6"


def test_event_trace_restart_truncates_to_checkpoint_and_continues_exactly(tmp_path):
    key = trace_entity_key("GB", "gb:1-2:0")
    recorder = EventTraceRecorder(
        tmp_path, run_id="run", pre_steps=1, post_steps=2, stride=1
    )
    recorder.capture(1, {key: _row("gb:1-2:0", 1.0)})
    recorder.trigger(
        {
            "event_id": "event-1",
            "event_type": "compatibility_release",
            "entity_id": "gb:1-2:0",
            "step": 2,
            "time": 2.0,
        },
        [key],
    )
    recorder.capture(2, {key: _row("gb:1-2:0", 2.0)})
    checkpoint = recorder.checkpoint_state()
    recorder.capture(3, {key: _row("gb:1-2:0", 999.0)})
    recorder.close()

    resumed = EventTraceRecorder(
        tmp_path,
        run_id="run",
        pre_steps=1,
        post_steps=2,
        stride=1,
        resume=True,
        checkpoint_state=checkpoint,
    )
    resumed.capture(3, {key: _row("gb:1-2:0", 3.0)})
    resumed.capture(4, {key: _row("gb:1-2:0", 4.0)})
    resumed.close()

    states = _read(tmp_path / "event_traces.csv")
    assert [float(row["local_normal_velocity"]) for row in states] == [1.0, 2.0, 3.0, 4.0]
    assert len(_read(tmp_path / "event_trace_events.csv")) == 1


def test_event_trace_restart_rejects_changed_window_shape(tmp_path):
    recorder = EventTraceRecorder(
        tmp_path, run_id="run", pre_steps=1, post_steps=2, stride=1
    )
    checkpoint = recorder.checkpoint_state()
    recorder.close()
    with pytest.raises(ValueError, match="window settings"):
        EventTraceRecorder(
            tmp_path,
            run_id="run",
            pre_steps=2,
            post_steps=2,
            stride=1,
            resume=True,
            checkpoint_state=checkpoint,
        )
