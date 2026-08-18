#!/usr/bin/env bash
set -euo pipefail

PROCESSES="${PROCESSES:-10}"
VIDEO_PROCESSES="${VIDEO_PROCESSES:-6}"
SKIP_TESTS="${SKIP_TESTS:-0}"
RUN_VIDEO="${RUN_VIDEO:-1}"
RUN_LONG="${RUN_LONG:-1}"
RENDER_VIDEO="${RENDER_VIDEO:-1}"

export PYTHONPATH="${PYTHONPATH:-src}"
mkdir -p results/production_summaries results/plots results/validation results/video_runs

python - <<'PY'
import importlib
mods = ["numpy", "numba", "scipy", "pyarrow", "yaml", "pandas", "matplotlib", "pytest"]
missing = []
for name in mods:
    try:
        importlib.import_module(name)
    except Exception as exc:
        missing.append(f"{name}: {exc}")
if missing:
    raise SystemExit("Missing runtime dependencies:\n  " + "\n  ".join(missing) +
                     "\nRun: python -m pip install -e '.[analysis,test]'")
print("Runtime dependency check passed.")
PY

echo "=== phase 1: baseline reconciliation + climb backpressure audit ==="
SKIP_TESTS="$SKIP_TESTS" PROCESSES="$PROCESSES" \
  bash scripts/run_next_audit_20260817.sh

if [[ "$RUN_VIDEO" == "1" ]]; then
  echo "=== phase 2: representative frame-preserving runs ==="
  VIDEO_ROOT=$(python scripts/run_video_cases.py \
    configs/production/video_representative_200.yaml \
    --processes "$VIDEO_PROCESSES" | tail -n 1)
  printf '%s\n' "$VIDEO_ROOT" | tee results/validation/overnight_video_root.txt

  if [[ "$RENDER_VIDEO" == "1" ]]; then
    echo "=== rendering representative videos ==="
    for run_dir in "$VIDEO_ROOT"/*; do
      [[ -d "$run_dir/frames" ]] || continue
      python scripts/render_microstructure_video.py "$run_dir" || true
    done
  fi
fi

if [[ "$RUN_LONG" == "1" ]]; then
  echo "=== phase 3: long-horizon selected-mechanism scaling ==="
  LONG_CAMPAIGN=$(python scripts/run_campaign.py \
    configs/production/overnight_long_horizon_200.yaml \
    --processes "$PROCESSES" --output-root results/campaigns | tail -n 1)
  printf '%s\n' "$LONG_CAMPAIGN" | tee results/validation/overnight_long_campaign_path.txt

  python scripts/analyze_campaign.py "$LONG_CAMPAIGN" \
    --output results/production_summaries/overnight_long_horizon_summary.csv
  python scripts/plot_campaign.py "$LONG_CAMPAIGN" \
    --summary results/production_summaries/overnight_long_horizon_summary.csv \
    --output results/plots/overnight_long_horizon
  python scripts/audit_barrierless_fraction.py "$LONG_CAMPAIGN" \
    --regime G2 --regime T2 --regime S2 --regime C5 \
    --output results/production_summaries/overnight_long_horizon_barrierless.csv
fi

echo "=== overnight campaign complete ==="
echo "Primary summaries: results/production_summaries/"
echo "Plots:             results/plots/"
echo "Video raw frames:  results/video_runs/"
