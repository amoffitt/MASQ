'''
Collapses tag to combine tags with one base error
Updates all counters to reflect rolled-up tags
'''
import argparse
import logging
import pickle
import sys

from typing import Optional

from masq.utils.io import tabprint
from masq.tags.cluster_rollup import cluster_rollup2


logger = logging.getLogger("masq_collapse_tags")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Collapses tag to combine tags with one base error"
            "Updates all counters to reflect rolled-up tags"
    )
    parser.add_argument(
        "--vt-counter",
        help="Pickle file containing the vt counter data for a single region",
    )
    parser.add_argument(
        "--vt-seq-counter",
        help="Pickle file containing the vt sequences counter data for a "
        "single region",
    )
    parser.add_argument(
        "--flip-counter",
        help="Pickle file containing the flip counter data for a "
        "single region",
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
        "--dna-input-ng",
        help="DNA input in nanograms",
        type=float,
    )

    parser.add_argument(
        "--output",
        help="Output tagcounts table filename",
    )
    parser.add_argument(
        "--output-vt-counter",
        help="Output vt counter pickle filename",
    )
    parser.add_argument(
        "--output-vt-seq-counter",
        help="Output vt seq counter pickle filename",
    )
    parser.add_argument(
        "--output-flip-counter",
        help="Output flip counter pickle filename",
    )
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    if argv is None:
        argv = sys.argv[1:]
    args = parse_args(argv)


    ########################################################################
    # Load sequence data
    print("loading sequence data...")
    with open(args.vt_counter, 'rb') as infile:
        vt_counter = pickle.load(infile)
    with open(args.vt_seq_counter, 'rb') as infile:
        vt_seq_counter = pickle.load(infile)
    with open(args.flip_counter, 'rb') as infile:
        flip_counter = pickle.load(infile)
    print("sequencing data loaded")

    ########################################################################
    # Do cluster rollup
    new_vt_counter, new_vt_seq_counter, _, _, _, new_flip_counter = \
        cluster_rollup2(
            vt_counter,
            vt_seq_counter,
            flip_counter,
            show_progress=False)

    ########################################################################

    # Save new counters
    with open(args.output_vt_counter, 'wb') as outfile:
        pickle.dump(new_vt_counter, outfile, pickle.HIGHEST_PROTOCOL)
    with open(args.output_vt_seq_counter, 'wb') as outfile:
        pickle.dump(new_vt_seq_counter, outfile, pickle.HIGHEST_PROTOCOL)
    with open(args.output_flip_counter, 'wb') as outfile:
        pickle.dump(new_flip_counter, outfile, pickle.HIGHEST_PROTOCOL)

    ########################################################################

    # Calculate tag counts for report
    num_obs_tags = sum(vt_counter.values())
    num_uniq_tags = len(vt_counter)
    num_collapsed_tags = len(new_vt_counter)

    ########################################################################
    # Report on original tags, unique tags, and rolled up tags
    with open(args.output, "w") as outfile:
        outfile.write(tabprint(["Region","Observed Tags", "Unique Tags", "Rolled-Up Tags",
                                "Avg Reads Per Unique Tag","Avg Reads Per Rolled-Up Tag",
                                "Fraction of Unique Tags that are Rolled-Up","Yield: Rolled-Up Tags"])+"\n")
        rolledupyield=float(num_collapsed_tags)/( args.dna_input_ng/(3.59*0.001) )
        counts = [int(args.region),num_obs_tags,num_uniq_tags,num_collapsed_tags,
                float(num_obs_tags)/max(1,num_uniq_tags), float(num_obs_tags)/max(1,num_collapsed_tags),
                float(num_collapsed_tags)/max(1,num_uniq_tags),rolledupyield]
        outfile.write(tabprint(counts)+"\n")
