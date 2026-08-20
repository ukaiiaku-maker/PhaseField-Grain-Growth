#!/usr/bin/env bash
set -euo pipefail

PROCESSES="${PROCESSES:-6}"
RUN_FULL_TESTS="${RUN_FULL_TESTS:-0}"

if [[ -n "${PYTHON_BIN:-}" ]]; then
  :
elif [[ -x /opt/anaconda3/bin/python ]]; then
  PYTHON_BIN=/opt/anaconda3/bin/python
elif command -v python3 >/dev/null 2>&1; then
  PYTHON_BIN="$(command -v python3)"
else
  PYTHON_BIN="$(command -v python)"
fi
export PYTHON_BIN
export PYTHONPATH="${PYTHONPATH:-src}"

mkdir -p results/preproduction_convergence results/production_summaries results/validation

echo "Using Python: $PYTHON_BIN"
"$PYTHON_BIN" - <<'PY'
import importlib, sys
for name in ("numpy", "numba", "scipy", "pyarrow", "yaml", "pandas", "matplotlib", "pytest"):
    importlib.import_module(name)
print("dependency check passed in", sys.executable)
PY

echo "=== physics/bookkeeping regression tests ==="
"$PYTHON_BIN" -m pytest -q \
  tests/integration/test_migration_closure.py \
  tests/unit/test_modes_climb_mechanics.py

if [[ "$RUN_FULL_TESTS" == "1" ]]; then
  echo "=== full regression suite ==="
  "$PYTHON_BIN" -m pytest -q
fi

# Re-analyze already-generated closure/video data if present. These consume no new
# simulation time and probe finite-population fitting and pixel-scale waviness.
CLOSURE_SOURCE="${CLOSURE_SOURCE:-}"
if [[ -z "$CLOSURE_SOURCE" ]]; then
  CLOSURE_SOURCE=$(ls -td results/migration_closure/* 2>/dev/null | head -n 1 || true)
fi
if [[ -n "$CLOSURE_SOURCE" && -f "$CLOSURE_SOURCE/campaign_manifest.json" ]]; then
  echo "=== existing-data population-window sensitivity ==="
  echo "$CLOSURE_SOURCE" | tee results/validation/preproduction_existing_closure_source.txt
  "$PYTHON_BIN" scripts/analyze_population_window_sensitivity.py "$CLOSURE_SOURCE" \
    --output results/production_summaries/preproduction_population_window_sensitivity.csv
else
  echo "No existing migration_closure campaign found; skipping population-window audit."
fi

VIDEO_SOURCE="${VIDEO_SOURCE:-}"
if [[ -z "$VIDEO_SOURCE" ]]; then
  VIDEO_SOURCE=$(ls -td results/migration_closure_video/* 2>/dev/null | head -n 1 || true)
fi
if [[ -n "$VIDEO_SOURCE" && -d "$VIDEO_SOURCE" ]]; then
  echo "=== existing-data boundary spectral roughness ==="
  echo "$VIDEO_SOURCE" | tee results/validation/preproduction_existing_video_source.txt
  "$PYTHON_BIN" scripts/analyze_boundary_spectral_roughness.py "$VIDEO_SOURCE" \
    --output results/production_summaries/preproduction_boundary_spectral_roughness.csv
else
  echo "No existing migration_closure_video campaign found; skipping spectral roughness audit."
fi

run_case () {
  local label="$1"
  local spec="$2"
  echo "=== convergence case: $label ===" >&2
  local campaign
  campaign=$("$PYTHON_BIN" scripts/run_migration_closure_campaign.py "$spec" \
    --processes "$PROCESSES" --output-root results/preproduction_convergence | tail -n 1)
  printf '%s\n' "$campaign" | tee "results/validation/preproduction_${label}_campaign.txt" >&2
  printf '%s\n' "$campaign"
}

REF=$(run_case reference configs/production/preproduction_convergence_reference_single.yaml | tail -n 1)
DT=$(run_case dt_half configs/production/preproduction_convergence_dt_half_single.yaml | tail -n 1)
GRID=$(run_case grid_fine configs/production/preproduction_convergence_grid_fine_single.yaml | tail -n 1)
SIZE=$(run_case size_large configs/production/preproduction_convergence_size_large_single.yaml | tail -n 1)

echo "=== actual-run gate-only bookkeeping audits ==="
for item in "reference:$REF" "dt_half:$DT" "grid_fine:$GRID" "size_large:$SIZE"; do
  label="${item%%:*}"
  path="${item#*:}"
  "$PYTHON_BIN" scripts/audit_gate_only_bookkeeping.py "$path" \
    --output "results/production_summaries/preproduction_bookkeeping_${label}.csv"
done

echo "=== normalized numerical-convergence analysis ==="
"$PYTHON_BIN" scripts/analyze_preproduction_convergence.py \
  --reference "$REF" \
  --dt-half "$DT" \
  --grid-fine "$GRID" \
  --size-large "$SIZE" \
  --output results/production_summaries/preproduction_convergence.csv

echo "=== preproduction freeze checks complete ==="
echo "reference campaign:  $REF"
echo "dt-half campaign:    $DT"
echo "fine-grid campaign:  $GRID"
echo "large-size campaign: $SIZE"
echo "convergence report:  results/production_summaries/preproduction_convergence.md"
echo "window report:       results/production_summaries/preproduction_population_window_sensitivity.md"
echo "roughness report:    results/production_summaries/preproduction_boundary_spectral_roughness.md"
