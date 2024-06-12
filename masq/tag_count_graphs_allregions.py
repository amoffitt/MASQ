'''
For all regions, plot number of reads per tag/template
'''
import argparse
import sys
from collections import Counter
import logging
import pickle
import time

from typing import Optional

import numpy as np

from masq.plots.tag_plots import (
    plot_number_of_reads_per_tag,
    plot_at_least_x_reads_per_tag,
)
from masq.utils.io import tabprint


logger = logging.getLogger("masq_tag_count_graphs_allgregions")

########################################################################

# Start timer
t0_all = time.time()


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="For all regions, plot number of reads per tag/template"
    )
    parser.add_argument(
        "--vt-counters",
        help="List of pickle files containing the vt counter data",
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

    ########################################################################
    # Load sequence data
    vt_counter_filenames = args.vt_counters.split(" ")
    vt_counters = []
    logger.debug("loading sequence data...")
    for vt_filename in vt_counter_filenames:
        with open(vt_filename, "rb") as infile:
            vt_counters.append(pickle.load(infile))
    logger.debug("sequencing data loaded")

    ########################################################################
    # tabulate distribution of reads per tag
    reads_per_tag: dict[str, int] = Counter()
    numreads_total = 0
    numtags_total = 0

    for vt_counter in vt_counters:
        tags = len(vt_counter)
        numtags_total += tags
        for _tag, tag_count in vt_counter.items():
            reads_per_tag[tag_count]+=1
            numreads_total += tag_count

    numreads = np.array(list(reads_per_tag.keys()))
    numtags = np.array(list(reads_per_tag.values()))

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
            outfile.write(tabprint(["Combined",reads,tags])+"\n")
