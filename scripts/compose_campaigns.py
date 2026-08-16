#!/usr/bin/env python3
from __future__ import annotations

import argparse

from grain_growth_pf.campaign import compose_completed_campaigns


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compose unique completed runs into one immutable analysis campaign."
    )
    parser.add_argument("campaigns", nargs="+")
    parser.add_argument("--output-root", default="results/campaigns")
    parser.add_argument("--expected-runs", type=int)
    arguments = parser.parse_args()
    print(compose_completed_campaigns(
        arguments.campaigns, arguments.output_root, arguments.expected_runs
    ))
