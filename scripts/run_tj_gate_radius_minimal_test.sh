#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-/opt/anaconda3/bin/python}"
PROCESSES="${PROCESSES:-1}"
REFERENCE_CAMPAIGN="${REFERENCE_CAMPAIGN:-results/preproduction_convergence/20260820T032540Z-8d9f46eb85}"
LEGACY_FINE_CAMPAIGN="${LEGACY_FINE_CAMPAIGN:-results/preproduction_convergence/20260820T074800Z-7180f14944}"
C5_VIDEO_ROOT="${C5_VIDEO_ROOT:-results/migration_closure_video/20260819T213538Z-76352e8475}"

export PYTHONPATH="${PYTHONPATH:-src}"
mkdir -p results/tj_gate_radius_minimal results/production_summaries results/validation

echo "Using Python: $PYTHON_BIN"
"$PYTHON_BIN" - <<'PY'
import importlib, sys
for name in ("numpy", "numba", "scipy", "pyarrow", "yaml", "pandas", "matplotlib"):
    importlib.import_module(name)
print("dependency check passed in", sys.executable)
PY

echo "=== re-render existing C5 movies with explicit mobility/pinning panel ==="
for name in C5_GATE_LINE-T900-s5101 C5_HYBRID_LINE-T900-s5101; do
  run="$C5_VIDEO_ROOT/$name"
  if [[ -d "$run/frames" ]]; then
    "$PYTHON_BIN" scripts/render_microstructure_video.py "$run" \
      --output "$run/pinning_footprint.gif" --no-png-frames --fps 20
  else
    echo "warning: no frames at $run; skipping"
  fi
done

echo "=== one new fine-grid T2 run with approximately fixed physical TJ footprint ==="
CORRECTED_ROOT=$("$PYTHON_BIN" scripts/run_migration_closure_video.py \
  configs/production/tj_gate_radius_minimal_fine_video.yaml \
  --output-root results/tj_gate_radius_minimal \
  --processes "$PROCESSES" | tail -n 1)
printf '%s\n' "$CORRECTED_ROOT" | tee results/validation/tj_gate_radius_minimal_video_root.txt

CORRECTED_RUN=$(find "$CORRECTED_ROOT" -maxdepth 1 -type d -name 'T2_GRID_FINE_TJ_R3-*' | head -n 1)
if [[ -z "$CORRECTED_RUN" ]]; then
  echo "could not locate corrected T2 run under $CORRECTED_ROOT" >&2
  exit 1
fi

echo "=== render corrected fine-grid T2 with blocked-domain and mobility panels ==="
"$PYTHON_BIN" scripts/render_microstructure_video.py "$CORRECTED_RUN" \
  --output "$CORRECTED_RUN/pinning_footprint.gif" --no-png-frames --fps 20

echo "=== focused kinetic comparison against already completed convergence runs ==="
"$PYTHON_BIN" scripts/analyze_tj_gate_radius_minimal.py \
  --reference-campaign "$REFERENCE_CAMPAIGN" \
  --legacy-fine-campaign "$LEGACY_FINE_CAMPAIGN" \
  --corrected-video-root "$CORRECTED_ROOT" \
  --output results/production_summaries/tj_gate_radius_minimal.csv

echo "=== minimal TJ gate test complete ==="
echo "corrected run: $CORRECTED_RUN"
echo "corrected movie: $CORRECTED_RUN/pinning_footprint.gif"
echo "existing gate-only C5 movie: $C5_VIDEO_ROOT/C5_GATE_LINE-T900-s5101/pinning_footprint.gif"
echo "existing hybrid C5 movie: $C5_VIDEO_ROOT/C5_HYBRID_LINE-T900-s5101/pinning_footprint.gif"
echo "analysis: results/production_summaries/tj_gate_radius_minimal.md"
