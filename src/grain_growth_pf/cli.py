from __future__ import annotations

import argparse

from grain_growth_pf.analysis.campaign import analyze_campaign
from grain_growth_pf.campaign import launch_campaign
from grain_growth_pf.config import ModelConfig
from grain_growth_pf.simulation import EventResolvedSimulation


def run_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    parser.add_argument("output")
    args = parser.parse_args()
    EventResolvedSimulation(ModelConfig.load(args.config), args.output).run()


def campaign_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("spec")
    parser.add_argument("--output-root", default="results/campaigns")
    parser.add_argument("--processes", type=int)
    args = parser.parse_args()
    print(launch_campaign(args.spec, args.output_root, args.processes))


def analyze_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--output")
    args = parser.parse_args()
    print(analyze_campaign(args.campaign, args.output).to_string(index=False))


def plot_main() -> None:
    from grain_growth_pf.analysis.plots import plot_campaign
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--output")
    args = parser.parse_args()
    print(plot_campaign(args.campaign, args.output))
