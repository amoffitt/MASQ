'''
Plot results of tag roll-up / collapsing for all loci
'''

import argparse
import logging
import sys
from typing import Optional

from masq.plots.rollup_results_plot import plot_rollup_results

logger = logging.getLogger("masq_plot_rollup_results")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot results of tag roll-up / collapsing for all loci"
    )

    parser.add_argument(
        "--combined-report",
        help="Rollup combined report filename",
    )

    parser.add_argument(
        "--plot1",
        help="Plot result of tag roll-up",
    )

    parser.add_argument(
        "--sample",
        help="Sample name/ID",
        default="",
    )

    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    if argv is None:
        argv = sys.argv[1:]

    args = parse_args(argv)

    input_file = args.combined_report
    output_file = args.plot1
    sample = args.sample
    ########################################################################

    # Process report file to get counts
    num_obs_tags = []
    num_uniq_tags = []
    num_collapsed_tags = []

    with open(input_file,"r") as infile:
        for linecount, line in enumerate(infile):
            parts = line.strip().split()
            if linecount==0:
                # skip the header
                continue

            _region = parts[0]
            num_obs_tags.append(int(parts[1]))
            num_uniq_tags.append(int(parts[2]))
            num_collapsed_tags.append(int(parts[3]))

    print(num_obs_tags)
    print(num_uniq_tags)
    print(num_collapsed_tags)

    ########################################################################
    # Graph for each region, total tags, unique tags, collapsed tags

    plot_rollup_results(
        sample,
        output_file,
        num_obs_tags,
        num_uniq_tags,
        num_collapsed_tags,
    )
