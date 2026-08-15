#!/usr/bin/env python3
import argparse

from grain_growth_pf.analysis.plots import plot_campaign


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--output")
    parser.add_argument("--summary")
    arguments = parser.parse_args()
    print(plot_campaign(arguments.campaign, arguments.output, arguments.summary))
