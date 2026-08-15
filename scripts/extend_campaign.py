#!/usr/bin/env python3
from __future__ import annotations

import argparse

from grain_growth_pf.campaign import extend_campaign


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Continue completed campaign runs from copied exact checkpoints."
    )
    parser.add_argument("source_campaigns", nargs="+", help="completed campaign directories")
    parser.add_argument("--max-steps", type=int, required=True)
    parser.add_argument("--termination-grains", type=int, required=True)
    parser.add_argument("--root", default="results/campaigns")
    parser.add_argument("--processes", type=int, default=None)
    args = parser.parse_args()
    result = extend_campaign(
        args.source_campaigns,
        args.max_steps,
        args.termination_grains,
        root=args.root,
        processes=args.processes,
    )
    print(result)


if __name__ == "__main__":
    main()
