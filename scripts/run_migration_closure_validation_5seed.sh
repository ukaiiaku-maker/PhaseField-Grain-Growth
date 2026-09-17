#!/usr/bin/env bash
set -euo pipefail

PROCESSES="${PROCESSES:-8}"
RUN_TESTS="${RUN_TESTS:-1}"

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
mkdir -p results/migration_closure_validation results/production_summaries results/plots results/validation

echo "Using Python: $PYTHON_BIN"
"$PYTHON_BIN" - <<'PY'
import importlib, sys
for name in ("numpy", "numba", "scipy", "pyarrow", "yaml", "pandas", "matplotlib", "pytest"):
    importlib.import_module(name)
print("dependency check passed in", sys.executable)
PY

if [[ "$RUN_TESTS" == "1" ]]; then
  echo "=== targeted migration-closure regression tests ==="
  "$PYTHON_BIN" -m pytest -q tests/integration/test_migration_closure.py
fi

echo "=== five-seed gate-only closure validation ==="
CAMPAIGN=$("$PYTHON_BIN" scripts/run_migration_closure_campaign.py \
  configs/production/migration_closure_validation_5seed.yaml \
  --processes "$PROCESSES" \
  --output-root results/migration_closure_validation | tail -n 1)
printf '%s\n' "$CAMPAIGN" | tee results/validation/migration_closure_validation_5seed_campaign.txt

STANDARD_SUMMARY=results/production_summaries/migration_closure_validation_5seed_standard.csv
"$PYTHON_BIN" scripts/analyze_campaign.py "$CAMPAIGN" --output "$STANDARD_SUMMARY"
"$PYTHON_BIN" scripts/plot_campaign.py "$CAMPAIGN" \
  --summary "$STANDARD_SUMMARY" \
  --output results/plots/migration_closure_validation_5seed

echo "=== matched-seed primary metrics ==="
"$PYTHON_BIN" scripts/analyze_migration_closure_validation.py "$CAMPAIGN" \
  --standard-summary "$STANDARD_SUMMARY" \
  --output-prefix results/production_summaries/migration_closure_validation_5seed

echo "=== effective-barrier audit ==="
"$PYTHON_BIN" scripts/audit_barrierless_fraction.py "$CAMPAIGN" \
  --regime G2_GATE_LINE \
  --regime T2_GATE_LINE \
  --regime S2_GATE_LINE \
  --regime C5_GATE_LINE \
  --regime C5_GATE_DIFFUSE \
  --output results/production_summaries/migration_closure_validation_5seed_barrierless.csv

echo "=== decision report ==="
cat results/production_summaries/migration_closure_validation_5seed_decision.md

echo "=== validation complete ==="
echo "campaign: $CAMPAIGN"
echo "standard summary: $STANDARD_SUMMARY"
echo "primary grouped metrics: results/production_summaries/migration_closure_validation_5seed_grouped.csv"
echo "paired seed metrics: results/production_summaries/migration_closure_validation_5seed_paired_to_B0.csv"
echo "decision report: results/production_summaries/migration_closure_validation_5seed_decision.md"
