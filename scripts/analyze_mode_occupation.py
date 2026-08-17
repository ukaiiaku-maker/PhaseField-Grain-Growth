#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from grain_growth_pf.analysis.mode_occupation import analyze_mode_occupation


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Measure temperature-dependent disconnection-mode occupation."
    )
    parser.add_argument("campaign")
    parser.add_argument("--regime", required=True)
    parser.add_argument("--bootstrap-samples", type=int, default=5000)
    parser.add_argument("--output", required=True)
    arguments = parser.parse_args()
    result = analyze_mode_occupation(
        arguments.campaign, arguments.regime, arguments.bootstrap_samples
    )
    target = Path(arguments.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2) + "\n")
    print(target)
