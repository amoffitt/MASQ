'''
Filters base tables to remove loci that have been filtered for QC metrics
QC fail loci defined in masq_QC_plots.R
'''

import argparse
import logging
import sys
from typing import Optional

logger = logging.getLogger("masq_qcfilter_report")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filters base tables to remove loci that have been "
        "filtered for QC metrics QC fail loci defined in masq_QC_plots.R"
    )

    parser.add_argument(
        "--qcfail-report",
        help="QC fail loci report filename",
    )
    parser.add_argument(
        "--base-report",
        help="Input combined base report filename",
    )
    parser.add_argument(
        "--filtered-base-report",
        help="Filtered base report filename to write in",
    )

    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    # Input is combined base report file
    # Output is same file, with QC filtered loci removed

    if argv is None:
        argv = sys.argv[1:]

    args = parse_args(argv)

    qc_fail_loci = set()
    with open(args.qcfail_report,'r') as f:
        for line in f:
            qc_fail_loci.add(line.strip())


    c=0
    with open(args.base_report,'r') as f:
        with open(args.filtered_base_report,'w') as fout:
            for line in f:
                if c == 0:
                    fout.write(line)
                else:
                    x=line.split()
                    if x[1] not in qc_fail_loci:
                        fout.write(line)
                c = c + 1
