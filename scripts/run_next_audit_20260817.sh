#!/usr/bin/env bash
set -euo pipefail

PROCESSES="${PROCESSES:-10}"
OUTPUT_ROOT="${OUTPUT_ROOT:-results/campaigns}"
SKIP_TESTS="${SKIP_TESTS:-0}"

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
mkdir -p results/production_summaries results/plots results/validation

echo "Using Python: $PYTHON_BIN"
"$PYTHON_BIN" - <<'PY'
import sys
import numpy, numba, scipy, pyarrow, yaml
print("python:", sys.executable)
print("numba:", numba.__version__)
PY

if [[ "$SKIP_TESTS" != "1" ]]; then
  echo "=== fast regression suite ==="
  "$PYTHON_BIN" -m pytest -q
fi

echo "=== baseline reconciliation: 20 matched seeds, B0/B1 ==="
BASE_CAMPAIGN=$("$PYTHON_BIN" scripts/run_campaign.py \
  configs/production/baseline_reconciliation_20.yaml \
  --processes "$PROCESSES" --output-root "$OUTPUT_ROOT" | tail -n 1)
printf '%s\n' "$BASE_CAMPAIGN" | tee results/validation/next_baseline_campaign_path.txt
"$PYTHON_BIN" scripts/analyze_campaign.py "$BASE_CAMPAIGN" \
  --output results/production_summaries/baseline_reconciliation_20_summary.csv
"$PYTHON_BIN" scripts/plot_campaign.py "$BASE_CAMPAIGN" \
  --summary results/production_summaries/baseline_reconciliation_20_summary.csv \
  --output results/plots/baseline_reconciliation_20

echo "=== climb/free-volume backpressure sweep ==="
CLIMB_CAMPAIGN=$("$PYTHON_BIN" scripts/run_campaign.py \
  configs/production/climb_backpressure_sweep_200.yaml \
  --processes "$PROCESSES" --output-root "$OUTPUT_ROOT" | tail -n 1)
printf '%s\n' "$CLIMB_CAMPAIGN" | tee results/validation/next_climb_campaign_path.txt
"$PYTHON_BIN" scripts/analyze_campaign.py "$CLIMB_CAMPAIGN" \
  --output results/production_summaries/climb_backpressure_sweep_200_summary.csv
"$PYTHON_BIN" scripts/plot_campaign.py "$CLIMB_CAMPAIGN" \
  --summary results/production_summaries/climb_backpressure_sweep_200_summary.csv \
  --output results/plots/climb_backpressure_sweep_200
"$PYTHON_BIN" scripts/audit_barrierless_fraction.py "$CLIMB_CAMPAIGN" \
  --output results/production_summaries/climb_backpressure_barrierless.csv

if [[ -n "${MECHANISM_CAMPAIGN:-}" ]]; then
  echo "=== attempt-limited audit of prior mechanism campaign ==="
  "$PYTHON_BIN" scripts/audit_barrierless_fraction.py "$MECHANISM_CAMPAIGN" \
    --regime E1 --regime SC3 --regime SC4 --regime J1 --regime J2 \
    --output results/production_summaries/mechanism_barrierless_audit.csv
fi

echo "=== complete ==="
echo "baseline campaign: $BASE_CAMPAIGN"
echo "climb campaign:    $CLIMB_CAMPAIGN"
echo "summaries: results/production_summaries/"
echo "plots:     results/plots/"
