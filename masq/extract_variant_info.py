'''
Extracts only the target variant bases from the full base report
'''

import argparse
import logging
import sys
from typing import Optional

logger = logging.getLogger("masq_extract_variant_info")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extracts only the target variant bases from the full "
        "base report"
    )

    parser.add_argument(
        "--combined-report",
        help="Input target variant bases report filename",
    )

    parser.add_argument(
        "--variant-report",
        help="Output variant report filename",
    )

    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    if argv is None:
        argv = sys.argv[1:]

    args = parse_args(argv)

    ########################################################################
    # Extract only the variant bases to a new report
    # Could add additional calculations here if necessary

    input_file = args.combined_report
    output_file = args.variant_report

    with open(input_file,"r") as infile, open(output_file,"w") as outfile:
        linecounter=1
        for x in infile:
            line = x.strip().split()
            if linecounter==1: # header
                outfile.write(x)
            elif int(line[0])== 2:
                outfile.write(x)
            linecounter+=1

    logger.info("Extracted variant bases to %s", output_file)
