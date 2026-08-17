#!/usr/bin/env python3
from __future__ import annotations

import argparse

from grain_growth_pf.analysis.pareto import write_pareto_summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Rank physically admissible jerkiness candidates by Pareto dominance."
    )
    parser.add_argument("summary")
    parser.add_argument("--output", required=True)
    parser.add_argument("--target-exponent", type=float, default=2.0)
    parser.add_argument("--exponent-tolerance", type=float, default=0.5)
    arguments = parser.parse_args()
    print(write_pareto_summary(
        arguments.summary, arguments.output,
        arguments.target_exponent, arguments.exponent_tolerance,
    ))
