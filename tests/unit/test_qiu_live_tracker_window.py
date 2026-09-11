import csv
import gzip
from pathlib import Path

from extract_qiu_live_tracker_window import extract


def _tracker(path: Path, rows: list[tuple[int, str]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["step", "value"])
        writer.writerows(rows)


def test_tracker_window_is_bounded_by_row_step(tmp_path):
    source = tmp_path / "live"
    source.mkdir()
    rows = [(9, "old"), (10, "first"), (11, "middle"), (12, "last"), (13, "new")]
    _tracker(source / "grain_tracks.csv", rows)
    _tracker(source / "boundary_tracks.csv", rows + [(14, "newer")])

    result = extract(source, tmp_path / "window", 10, 12)

    assert result["files"]["grain_tracks.csv"]["rows"] == 3
    assert result["files"]["boundary_tracks.csv"]["step_last"] == 12
    with gzip.open(tmp_path / "window" / "grain_tracks.csv.gz", "rt", newline="") as stream:
        assert [int(row["step"]) for row in csv.DictReader(stream)] == [10, 11, 12]
