#!/usr/bin/env bash
set -euo pipefail

PROCESSES="${PROCESSES:-10}"
VIDEO_PROCESSES="${VIDEO_PROCESSES:-6}"
SKIP_TESTS="${SKIP_TESTS:-0}"
RUN_VIDEO="${RUN_VIDEO:-1}"
RUN_LONG="${RUN_LONG:-1}"
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

export PYTHON_BIN
export PYTHONPATH="${PYTHONPATH:-src}"
mkdir -p results/production_summaries results/plots results/validation results/video_runs

echo "Using Python: $PYTHON_BIN"
"$PYTHON_BIN" - <<'PY'
import importlib, sys
mods = ["numpy", "numba", "scipy", "pyarrow", "yaml", "pandas", "matplotlib", "pytest"]
missing = []
for name in mods:
    try:
        importlib.import_module(name)
    except Exception as exc:
        missing.append(f"{name}: {exc}")
if missing:
    raise SystemExit("Missing runtime dependencies in " + sys.executable + ":\n  " + "\n  ".join(missing))
print("Runtime dependency check passed in", sys.executable)
PY

echo "=== phase 1: baseline reconciliation + climb backpressure audit ==="
PYTHON_BIN="$PYTHON_BIN" SKIP_TESTS="$SKIP_TESTS" PROCESSES="$PROCESSES" \
  bash scripts/run_next_audit_20260817.sh

if [[ "$RUN_VIDEO" == "1" ]]; then
  echo "=== phase 2: representative frame-preserving runs ==="
  VIDEO_ROOT=$("$PYTHON_BIN" scripts/run_video_cases.py \
    configs/production/video_representative_200.yaml \
    --processes "$VIDEO_PROCESSES" | tail -n 1)
  printf '%s\n' "$VIDEO_ROOT" | tee results/validation/overnight_video_root.txt

  if [[ "$RENDER_VIDEO" == "1" ]]; then
    echo "=== rendering representative PNG sequences and videos ==="
    "$PYTHON_BIN" scripts/render_video_campaign.py "$VIDEO_ROOT" || true
  fi
fi

if [[ "$RUN_LONG" == "1" ]]; then
  echo "=== phase 3: long-horizon selected-mechanism scaling ==="
  LONG_CAMPAIGN=$("$PYTHON_BIN" scripts/run_campaign.py \
    configs/production/overnight_long_horizon_200.yaml \
    --processes "$PROCESSES" --output-root results/campaigns | tail -n 1)
  printf '%s\n' "$LONG_CAMPAIGN" | tee results/validation/overnight_long_campaign_path.txt

  "$PYTHON_BIN" scripts/analyze_campaign.py "$LONG_CAMPAIGN" \
    --output results/production_summaries/overnight_long_horizon_summary.csv
  "$PYTHON_BIN" scripts/plot_campaign.py "$LONG_CAMPAIGN" \
    --summary results/production_summaries/overnight_long_horizon_summary.csv \
    --output results/plots/overnight_long_horizon
  "$PYTHON_BIN" scripts/audit_barrierless_fraction.py "$LONG_CAMPAIGN" \
    --regime G2 --regime T2 --regime S2 --regime C5 \
    --output results/production_summaries/overnight_long_horizon_barrierless.csv
fi

echo "=== overnight campaign complete ==="
echo "Primary summaries: results/production_summaries/"
echo "Plots:             results/plots/"
echo "Video outputs:     results/video_runs/"
