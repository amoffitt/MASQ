'''
Combines within-template error statistics across all loci
'''
import argparse
import sys
from typing import Optional
import pickle

import numpy as np

from masq.utils.io import tabprint

def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combines within-template error statistics across "
        "all loci"
    )
    parser.add_argument(
        "--within-tag-errors",
        help="Pickle files for within tag errors",
    )
    parser.add_argument(
        "--output-table1",
        help="Output table filename",
    )
    parser.add_argument(
        "--output-table2",
        help="Output table filename",
    )

    return parser.parse_args(argv)

########################################################################
NUCS = ["A", "C", "G", "T"]
REF_TRINUCS = [x+y+z for x in NUCS for y in NUCS for z in NUCS]
PAIRED_TRINUCS=[ (a , a[0]+x+a[2]) for a in REF_TRINUCS for x in NUCS if x!=a[1] ]
ERR_RANGE1 = np.arange(-0.1,1.0,0.1)

# Intiialize table


def main(argv: Optional[list[str]] = None) -> None:
    if argv is None:
        argv = sys.argv[1:]

    args = parse_args(argv)
    input_file_list = args.within_tag_errors.split(" ")

    within_tag_errors_all = np.zeros(
        shape=(len(PAIRED_TRINUCS),len(ERR_RANGE1),2), dtype=int)

    for f in input_file_list:
        # Add to current table
        with open(f, "rb") as infile:
            within_tag_errs_region = pickle.load(infile)
        within_tag_errors_all = within_tag_errors_all  +  within_tag_errs_region

    with open(args.output_table1, "w") as outfile1, \
            open(args.output_table2, "w") as outfile2:
        # Original table
        for i,pt in enumerate(PAIRED_TRINUCS):
            outfile1.write(tabprint(
                ['R1']+list(pt)+list(within_tag_errors_all[i,:,0]))+'\n')
            outfile1.write(tabprint(
                ['R2']+list(pt)+list(within_tag_errors_all[i,:,1]))+'\n')

        # Convert to fractions, skip the 0 error bin
        for i,pt in enumerate(PAIRED_TRINUCS):
            frac_r1 = within_tag_errors_all[i,1:,0] / \
                within_tag_errors_all[i,1:,0].sum(keepdims=True)
            frac_r2 = within_tag_errors_all[i,1:,1] / \
                within_tag_errors_all[i,1:,1].sum(keepdims=True)
            outfile2.write(tabprint(['R1']+list(pt)+list(frac_r1))+'\n')
            outfile2.write(tabprint(['R2']+list(pt)+list(frac_r2))+'\n')
