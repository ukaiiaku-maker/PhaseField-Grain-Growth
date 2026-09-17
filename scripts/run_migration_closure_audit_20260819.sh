#!/usr/bin/env bash
set -euo pipefail

PROCESSES="${PROCESSES:-8}"
RUN_FULL_TESTS="${RUN_FULL_TESTS:-1}"
RENDER_VIDEO="${RENDER_VIDEO:-1}"

if [[ -n "${PYTHON_BIN:-}" ]]; then
  :
elif [[ -x /opt/anaconda3/bin/python ]]; then
  PYTHON_BIN=/opt/anaconda3/bin/python
elif command -v python3 >/dev/null 2>&1; then
  PYTHON_BIN="$(command -v python3)"
else
  PYTHON_BIN="$(command -v python)"
fi
export PYTHONPATH="${PYTHONPATH:-src}"

mkdir -p results/validation results/production_summaries results/plots

echo "Using Python: $PYTHON_BIN"
"$PYTHON_BIN" - <<'PY'
import numpy, numba, scipy, pyarrow, pandas, matplotlib, pytest
print("dependency check passed")
PY

echo "=== targeted migration-closure tests ==="
"$PYTHON_BIN" -m pytest -q tests/integration/test_migration_closure.py

if [[ "$RUN_FULL_TESTS" == "1" ]]; then
  echo "=== full fast regression suite ==="
  "$PYTHON_BIN" -m pytest -q
fi

echo "=== matched closure smoke campaign ==="
CLOSURE_ROOT=$("$PYTHON_BIN" scripts/run_migration_closure_campaign.py \
  configs/production/migration_closure_smoke_200.yaml \
  --processes "$PROCESSES" | tail -n 1)
printf '%s\n' "$CLOSURE_ROOT" | tee results/validation/migration_closure_smoke_path.txt

"$PYTHON_BIN" scripts/analyze_campaign.py "$CLOSURE_ROOT" \
  --output results/production_summaries/migration_closure_smoke_summary.csv
"$PYTHON_BIN" scripts/plot_campaign.py "$CLOSURE_ROOT" \
  --summary results/production_summaries/migration_closure_smoke_summary.csv \
  --output results/plots/migration_closure_smoke

echo "=== matched C5 geometry/video comparison ==="
VIDEO_ROOT=$("$PYTHON_BIN" scripts/run_migration_closure_video.py \
  configs/production/migration_closure_video_c5.yaml \
  --processes 4 | tail -n 1)
printf '%s\n' "$VIDEO_ROOT" | tee results/validation/migration_closure_video_path.txt

"$PYTHON_BIN" scripts/analyze_frame_roughness.py "$VIDEO_ROOT" \
  --output results/production_summaries/migration_closure_boundary_roughness.csv

if [[ "$RENDER_VIDEO" == "1" ]]; then
  for run_dir in "$VIDEO_ROOT"/*; do
    [[ -d "$run_dir/frames" ]] || continue
    "$PYTHON_BIN" scripts/render_microstructure_video.py "$run_dir" || true
  done
fi

echo "=== migration-closure audit complete ==="
echo "closure campaign: $CLOSURE_ROOT"
echo "video campaign:   $VIDEO_ROOT"
echo "summary: results/production_summaries/migration_closure_smoke_summary.csv"
echo "roughness: results/production_summaries/migration_closure_boundary_roughness.csv"
