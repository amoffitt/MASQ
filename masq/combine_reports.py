'''
Combines reports from individual loci into one report
'''
import argparse
import logging
import sys
from typing import Optional

logger = logging.getLogger("masq_combine_reports")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combines reports from individual loci into one report"
    )

    parser.add_argument(
        "--region-reports",
        help="List of input region reports filenames",
    )

    parser.add_argument(
        "--combined-report",
        help="Output combined report filename",
    )

    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    if argv is None:
        argv = sys.argv[1:]

    args = parse_args(argv)

    input_file_list = args.region_reports.split(" ")
    output_file = args.combined_report

    with open(output_file,"w") as outfile:

        counter=0
        for input_filename in input_file_list:
            counter +=1
            with open(input_filename, "r") as infile:
                linecounter=1
                for line in infile:
                    if linecounter==1:
                        if counter==1:
                            # first file header
                            outfile.write(line)
                    else:
                        outfile.write(line)
                    linecounter+=1
