#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-/opt/anaconda3/bin/python}"
PROCESSES="${PROCESSES:-4}"
RENDER="${RENDER:-1}"

export PYTHONPATH="${PYTHONPATH:-src}"
mkdir -p results/jerky_integrity_smoke results/production_summaries results/validation

echo "Using Python: $PYTHON_BIN"
"$PYTHON_BIN" - <<'PY'
import importlib, sys
for name in ("numpy", "numba", "scipy", "pyarrow", "yaml", "pandas", "matplotlib", "pytest"):
    importlib.import_module(name)
print("dependency check passed in", sys.executable)
PY

echo "=== focused jerky-growth regression tests ==="
"$PYTHON_BIN" -m pytest -q \
  tests/geometry/test_arclength_tracker.py \
  tests/integration/test_migration_closure.py \
  tests/integration/test_jerky_mechanism_integrity.py \
  tests/unit/test_activation_work_decomposition.py \
  tests/unit/test_modes_climb_mechanics.py

echo "=== short B0/G/T/GT/S/C/coupled/Qiu integrity campaign ==="
ROOT=$("$PYTHON_BIN" scripts/run_migration_closure_video.py \
  configs/production/jerky_mechanism_integrity_smoke.yaml \
  --output-root results/jerky_integrity_smoke \
  --processes "$PROCESSES" | tail -n 1)
printf '%s\n' "$ROOT" | tee results/validation/jerky_integrity_smoke_root.txt

echo "=== scaling, jerkiness, shear, and activation-work audit ==="
"$PYTHON_BIN" scripts/analyze_jerky_mechanism_integrity.py "$ROOT" \
  --output results/production_summaries/jerky_mechanism_integrity_smoke.csv

if [[ "$RENDER" == "1" ]]; then
  echo "=== render selected mechanism movies ==="
  for regime in B0 G T GT S GS TS GTS C GTC GTSC QIU; do
    run=$(find "$ROOT" -maxdepth 1 -type d -name "${regime}-T900-s5101" | head -n 1)
    if [[ -n "$run" ]]; then
      "$PYTHON_BIN" scripts/render_microstructure_video.py "$run" \
        --output "$run/integrity.gif" --no-png-frames --fps 20
    fi
  done
fi

echo "=== jerky mechanism smoke complete ==="
echo "root: $ROOT"
echo "analysis: results/production_summaries/jerky_mechanism_integrity_smoke.md"
