#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-/opt/anaconda3/bin/python}"
PROCESSES="${PROCESSES:-6}"

export PYTHONPATH="${PYTHONPATH:-src}"
mkdir -p results/jerky_discriminating_900K results/production_summaries results/validation

echo "Using Python: $PYTHON_BIN"
echo "=== G/T/GT barrier, GS packet, and C activation-scale screen at 900 K ==="
ROOT=$("$PYTHON_BIN" scripts/run_migration_closure_campaign.py \
  configs/production/jerky_mechanism_discriminating_900K.yaml \
  --output-root results/jerky_discriminating_900K \
  --processes "$PROCESSES" | tail -n 1)
printf '%s\n' "$ROOT" | tee results/validation/jerky_discriminating_900K_root.txt

"$PYTHON_BIN" scripts/analyze_jerky_mechanism_integrity.py "$ROOT" \
  --output results/production_summaries/jerky_mechanism_discriminating_900K.csv

echo "=== discriminating mechanism screen complete ==="
echo "root: $ROOT"
echo "analysis: results/production_summaries/jerky_mechanism_discriminating_900K.md"
