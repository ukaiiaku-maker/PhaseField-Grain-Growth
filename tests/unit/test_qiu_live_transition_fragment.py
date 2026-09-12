from argparse import Namespace
import json
from pathlib import Path
import tarfile

import pyarrow as pa
import pyarrow.parquet as pq

from create_qiu_live_transition_fragment import create_fragment


def test_transition_fragment_retains_bounded_closed_rows_and_dense_window(tmp_path):
    source = tmp_path / "legacy-run"
    source.mkdir()
    checkpoint = {
        "step_number": 10,
        "time": 1.0,
        "grain_tracks_offset": 3,
        "boundary_tracks_offset": 4,
    }
    (source / "checkpoint.json").write_text(json.dumps(checkpoint))
    (source / "checkpoint.npz").write_bytes(b"checkpoint")
    (source / "manifest.json").write_text("{}")
    (source / "grain_tracks.csv").write_bytes(b"abcdef")
    (source / "boundary_tracks.csv").write_bytes(b"abcdefgh")
    parts = source / "per_step_diagnostics.parquet"
    parts.mkdir()
    for index, steps in enumerate(((11, 12), (13, 14), (15,), (16,))):
        pq.write_table(
            pa.table({"step": pa.array(steps, type=pa.int64())}),
            parts / f"part-{index:06d}.parquet",
        )
    boundary_parts = source / "per_boundary_diagnostics.parquet"
    boundary_parts.mkdir()
    for index, steps in enumerate(((11, 12, 13), (14, 15), (16,))):
        pq.write_table(
            pa.table({"step": pa.array(steps, type=pa.int64())}),
            boundary_parts / f"part-{index:06d}.parquet",
        )
    fields = source / "diagnostic_fields"
    fields.mkdir()
    for step in range(11, 17):
        (fields / f"step-{step:07d}-guard.npz").write_bytes(str(step).encode())
    frames = source / "frames"
    frames.mkdir()
    (frames / "frame-0000010.npz").write_bytes(b"10")
    (frames / "frame-0000020.npz").write_bytes(b"20")
    (frames / "._frame-0000010.npz").write_bytes(b"appledouble")

    result = create_fragment(Namespace(
        source=source,
        destination=tmp_path / "output",
        archive_name="fragment.tar.gz",
        checkpoint_step=10,
        through_step=15,
        closed_part_count=3,
        closed_boundary_part_count=2,
        diagnostic_start=11,
        coarse_field_stride=2,
        dense_field_start=14,
        scratch=tmp_path,
    ))

    assert result["diagnostics_complete_through_step"] == 15
    assert result["closed_scalar_parts"] == 3
    assert result["closed_boundary_parts"] == 2
    assert result["scalar_stream"] == {
        "row_count": 5,
        "unique_step_count": 5,
        "step_first": 11,
        "step_last": 15,
    }
    assert result["boundary_stream"]["step_last"] == 15
    assert result["field_unique_step_count"] == 4
    with tarfile.open(tmp_path / "output" / "fragment.tar.gz") as archive:
        names = set(archive.getnames())
        prefix = "qualification/legacy-run"
        assert f"{prefix}/per_step_diagnostics.parquet/part-000002.parquet" in names
        assert f"{prefix}/per_step_diagnostics.parquet/part-000003.parquet" not in names
        assert f"{prefix}/per_boundary_diagnostics.parquet/part-000001.parquet" in names
        assert f"{prefix}/per_boundary_diagnostics.parquet/part-000002.parquet" not in names
        assert f"{prefix}/diagnostic_fields/step-0000011-guard.npz" in names
        assert f"{prefix}/diagnostic_fields/step-0000012-guard.npz" in names
        assert f"{prefix}/diagnostic_fields/step-0000013-guard.npz" not in names
        assert f"{prefix}/diagnostic_fields/step-0000014-guard.npz" in names
        assert f"{prefix}/diagnostic_fields/step-0000015-guard.npz" in names
        assert f"{prefix}/frames/frame-0000010.npz" in names
        assert f"{prefix}/frames/._frame-0000010.npz" not in names
        assert f"{prefix}/frames/frame-0000020.npz" not in names


def test_transition_fragment_rejects_incorrect_assumed_part_boundary(tmp_path):
    source = tmp_path / "legacy-run"
    source.mkdir()
    checkpoint = {
        "step_number": 10,
        "time": 1.0,
        "grain_tracks_offset": 0,
        "boundary_tracks_offset": 0,
    }
    (source / "checkpoint.json").write_text(json.dumps(checkpoint))
    (source / "checkpoint.npz").write_bytes(b"checkpoint")
    (source / "manifest.json").write_text("{}")
    (source / "grain_tracks.csv").write_bytes(b"")
    (source / "boundary_tracks.csv").write_bytes(b"")
    parts = source / "per_step_diagnostics.parquet"
    parts.mkdir()
    pq.write_table(pa.table({"step": [11, 12]}), parts / "part-000000.parquet")
    (source / "per_boundary_diagnostics.parquet").mkdir()
    (source / "diagnostic_fields").mkdir()
    (source / "frames").mkdir()

    import pytest

    with pytest.raises(RuntimeError, match="end at step 12"):
        create_fragment(Namespace(
            source=source,
            destination=tmp_path / "output",
            archive_name="fragment.tar.gz",
            checkpoint_step=10,
            through_step=13,
            closed_part_count=1,
            closed_boundary_part_count=0,
            diagnostic_start=11,
            coarse_field_stride=2,
            dense_field_start=12,
            scratch=tmp_path,
        ))
