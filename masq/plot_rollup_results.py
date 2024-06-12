'''
Plot results of tag roll-up / collapsing for all loci
'''

import os
import numpy as np
import gzip
from collections import Counter, defaultdict
import fileinput
import operator
import pickle
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm


import argparse
import logging
import sys
from typing import Optional

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


WIDTH = 0.25


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
        linecount=1
        regioncount = 0
        for line in infile:
            parts = line.strip().split()
            if linecount>1:
                region = parts[0]
                regioncount+=1
                num_obs_tags.append(int(parts[1]))
                num_uniq_tags.append(int(parts[2]))
                num_collapsed_tags.append(int(parts[3]))
            linecount+=1

    print(num_obs_tags)
    print(num_uniq_tags)
    print(num_collapsed_tags)

    ########################################################################

    # Graph for each region, total tags, unique tags, collapsed tags
    N = regioncount
    fig = plt.figure(figsize=(50,10))
    fig.suptitle(sample+"\n"+
                "Results of Tag Rollup - 1 error allowed", fontsize=16)
    x = list(range(regioncount))
    ax = plt.subplot(111)
    ax.bar(
        x, num_obs_tags, WIDTH,
        alpha=0.5, color='purple', label='Total')
    ax.bar(
        [p + WIDTH for p in x], num_uniq_tags, WIDTH,
        alpha=0.5, color='green', label='Unique')
    ax.bar(
        [p + WIDTH*2 for p in x], num_collapsed_tags, WIDTH,
        alpha=0.5, color='blue', label='Collapsed')
    ax.legend(['Total','Unique','Collapsed'], loc='upper left')
    ax.ticklabel_format(style='plain')
    ax.get_yaxis().set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

    plt.xticks([p + WIDTH for p in x],range(N))
    plt.xlim([min(x)-WIDTH,max(x)+WIDTH*4])
    plt.xlabel("Region",fontsize=14)
    plt.ylabel("Tag Counts",fontsize=14)
    plt.savefig(output_file, dpi=200, facecolor='w', edgecolor='w',
                papertype=None, format=None,
                transparent=False)
    plt.close()
