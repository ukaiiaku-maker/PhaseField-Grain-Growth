from argparse import Namespace
from pathlib import Path
import json
import tarfile

from create_qiu_live_snapshot import create_snapshot


def _write_parts(root: Path, count: int) -> None:
    root.mkdir(parents=True)
    for index in range(count):
        (root / f"part-{index:06d}.parquet").write_bytes(str(index).encode())


def test_live_snapshot_uses_checkpoint_counters_and_offsets(tmp_path: Path) -> None:
    source = tmp_path / "qualification"
    source.mkdir()
    checkpoint = {
        "step_number": 100,
        "time": 4.0,
        "grain_tracks_offset": 3,
        "boundary_tracks_offset": 4,
        "extension_state": {
            "qiu_forensics": {
                "scalar_part": 2,
                "boundary_part": 1,
                "population_history": [[100, 42]],
            },
            "event_trace": {"state_part_count": 1, "event_part_count": 1},
        },
    }
    (source / "checkpoint.json").write_text(json.dumps(checkpoint))
    (source / "checkpoint.npz").write_bytes(b"checkpoint")
    (source / "manifest.json").write_text("{}")
    (source / "grain_tracks.csv").write_bytes(b"abcdef")
    (source / "boundary_tracks.csv").write_bytes(b"abcdefgh")
    _write_parts(source / "per_step_diagnostics.parquet", 3)
    _write_parts(source / "per_boundary_diagnostics.parquet", 2)
    _write_parts(source / "event_traces.parquet", 2)
    _write_parts(source / "event_trace_events.parquet", 2)
    (source / "events.parquet").mkdir()
    (source / "frames").mkdir()
    (source / "frames" / "frame-0000000.npz").write_bytes(b"0")
    (source / "frames" / "frame-0000100.npz").write_bytes(b"100")
    (source / "frames" / "frame-0000200.npz").write_bytes(b"200")
    (source / "diagnostic_fields").mkdir()
    (source / "diagnostic_fields" / "step-0000100-cadence.npz").write_bytes(b"100")
    (source / "diagnostic_fields" / "step-0000101-cadence.npz").write_bytes(b"101")

    result = create_snapshot(
        Namespace(
            source=source,
            destination=tmp_path / "durable",
            archive_name="snapshot.tar.gz",
            expected_step=100,
            layout="corrected",
            scratch=tmp_path,
        )
    )

    assert result["closed_scalar_parts"] == 2
    assert result["last_grain_count"] == 42
    with tarfile.open(tmp_path / "durable" / "snapshot.tar.gz") as archive:
        names = set(archive.getnames())
        assert "qualification/per_step_diagnostics.parquet/part-000001.parquet" in names
        assert "qualification/per_step_diagnostics.parquet/part-000002.parquet" not in names
        assert "qualification/frames/frame-0000100.npz" in names
        assert "qualification/frames/frame-0000200.npz" not in names
        assert archive.extractfile("qualification/grain_tracks.csv").read() == b"abc"
        assert archive.extractfile("qualification/boundary_tracks.csv").read() == b"abcd"

