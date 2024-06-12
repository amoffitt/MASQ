'''
For a single region, plot the number of reads per template/tag
'''
import argparse
import logging
from collections import Counter
import pickle
from typing import Optional
import sys

import numpy as np

from masq.utils.io import tabprint

from masq.plots.tag_plots import (
    plot_number_of_reads_per_tag,
    plot_at_least_x_reads_per_tag,
)

logger = logging.getLogger("masq_tag_count_graphs")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="For a single region, plot the number of reads per "
        "template/tag"
    )
    parser.add_argument(
        "--vt-counter",
        help="Pickle file containing the vt counter data for a single region",
    )
    parser.add_argument(
        "--region",
        help="Region number",
    )
    parser.add_argument(
        "--sample",
        help="Sample name/ID",
        default="",
    )
    parser.add_argument(
        "--output",
        help="Output tagcounts table filename",
    )

    parser.add_argument(
        "--plot1",
        help="Number of reads per tag plot filename",
    )
    parser.add_argument(
        "--plot2",
        help="Atleast reads per tag plot filename",
    )

    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    if argv is None:
        argv = sys.argv[1:]

    args = parse_args(argv)

    logger.info("loading sequence data...")
    with open(args.vt_counter, "rb") as infile:
        vt_counter = pickle.load(infile)
    logger.info("sequence data loaded")

    ########################################################################
    # tabulate distribution of reads per tag
    reads_per_tag: dict[str, int] = Counter()
    numreads_total = 0
    numtags_total = len(vt_counter)
    for _tag, tag_count in vt_counter.items():
        reads_per_tag[tag_count]+=1
        numreads_total += tag_count
    numreads=np.array(list(reads_per_tag.keys()))
    numtags=np.array(list(reads_per_tag.values()))

    ########################################################################
    # Plot reads per tag
    plot_number_of_reads_per_tag(
        numreads,
        numtags,
        numreads_total,
        numtags_total,
        sample=args.sample,
        filename=args.plot1,
        logscale=True)
    plot_at_least_x_reads_per_tag(
        numreads,
        numtags,
        numreads_total,
        numtags_total,
        sample=args.sample,
        filename=args.plot2)

    ########################################################################
    # Write bar graph counts out to file
    with open(args.output, "w") as outfile:
        header = ["Region","Reads Per Tag","Number of Tags"]
        outfile.write(tabprint(header)+"\n")
        for reads, tags in zip(numreads, numtags):
            outfile.write(tabprint([int(args.region),reads,tags])+"\n")
