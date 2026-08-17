#!/usr/bin/env python3
from __future__ import annotations

import argparse

from grain_growth_pf.analysis.aggregate import aggregate_summaries


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Aggregate finalized campaign summaries; later inputs take precedence."
    )
    parser.add_argument("summaries", nargs="+")
    parser.add_argument("--output", required=True)
    arguments = parser.parse_args()
    print(aggregate_summaries(arguments.summaries, arguments.output))
