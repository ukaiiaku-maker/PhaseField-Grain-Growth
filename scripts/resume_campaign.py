#!/usr/bin/env python3
from __future__ import annotations

import argparse

from grain_growth_pf.campaign import resume_campaign


def main() -> None:
    parser = argparse.ArgumentParser(description="Resume an interrupted campaign in place")
    parser.add_argument("campaign")
    parser.add_argument("--processes", type=int)
    args = parser.parse_args()
    print(resume_campaign(args.campaign, args.processes))


if __name__ == "__main__":
    main()
