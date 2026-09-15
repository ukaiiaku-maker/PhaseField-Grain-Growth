#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-/opt/anaconda3/bin/python}"
PROCESSES="${PROCESSES:-6}"

export PYTHONPATH="${PYTHONPATH:-src}"
mkdir -p results/jerky_factorial_900K results/production_summaries results/validation

echo "Using Python: $PYTHON_BIN"

echo "=== one-seed full B0/G/T/S/C factorial campaign at 900 K ==="
ROOT=$("$PYTHON_BIN" scripts/run_migration_closure_campaign.py \
  configs/production/jerky_mechanism_factorial_900K.yaml \
  --output-root results/jerky_factorial_900K \
  --processes "$PROCESSES" | tail -n 1)
printf '%s\n' "$ROOT" | tee results/validation/jerky_factorial_900K_root.txt

"$PYTHON_BIN" scripts/analyze_jerky_mechanism_integrity.py "$ROOT" \
  --output results/production_summaries/jerky_mechanism_factorial_900K.csv

echo "=== factorial campaign complete ==="
echo "root: $ROOT"
echo "analysis: results/production_summaries/jerky_mechanism_factorial_900K.md"
