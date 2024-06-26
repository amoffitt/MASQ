'''
Combines reports from different steps of the pipeline into one final report per sample
'''

import argparse
import logging
import sys
from typing import Optional

from masq.utils.io import tabprint

logger = logging.getLogger("masq_final_report")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combines reports from different steps of the pipeline "
        "into one final report per sample"
    )

    parser.add_argument(
        "--input-snv-table",
        help="Extended SNV table filename",
    )
    parser.add_argument(
        "--report-primers",
        help="Primer counters report filename",
    )
    parser.add_argument(
        "--report-rollup",
        help="Rollup report filename",
    )
    parser.add_argument(
        "--report-alignment",
        help="Alignment counters report filename",
    )
    parser.add_argument(
        "--report-variants",
        help="Variant base counters report filename",
    )
    parser.add_argument(
        "--combined-report",
        help="Output combined report filename",
    )

    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    """Combines reports from different steps of the pipeline into one report.

    Combines reports from different steps of the pipeline into one final
    report per sample
    """
    if argv is None:
        argv = sys.argv[1:]

    args = parse_args(argv)

    # Input is list of region specific report files which all have a header
    # Output is combined report files with one header
    logger.info('Combining reports')
    header_all = []
    data_all = []

    with open(args.input_snv_table,'r') as f:
        header = f.readline().strip("\n").split("\t")
        header_all.extend(header)
        for line in f:
            x = line.strip("\n").split("\t")
            data_all.append(x)

    print(data_all)

    with open(args.report_primers,'r') as f:
        c=0
        _header1 = f.readline()
        _summarydata = f.readline()
        f.readline()
        header=f.readline().strip("\n").split("\t")[1:]
        header_all.extend(header)
        for line in f:
            x = line.strip("\n").split("\t")[1:]
            if not line.startswith("NONE"):
                data_all[c].extend(x)
                c+=1

    with open(args.report_rollup,'r') as f:
        c=0
        header = f.readline().strip("\n").split("\t")[1:]
        header_all.extend(header)
        for line in f:
            x = line.strip("\n").split("\t")[1:]
            data_all[c].extend(x)
            c+=1

    with open(args.report_alignment,'r') as f:
        header = f.readline().strip("\n").split("\t")[2:]
        header_all.extend(header)
        align_data: dict = {}
        regct=0
        for line in f:
            x = line.strip("\n").split("\t")
            reg=int(x[0])
            if reg in align_data:
                align_data[reg] = [int(x[2]),int(x[3]),float(x[4]),int(x[5]),int(x[6]) + align_data[reg][4] ,float(x[7]) + align_data[reg][5]]
            else:
                align_data[reg] = [int(x[2]),int(x[3]),float(x[4]),int(x[5]),int(x[6]),float(x[7])]
                regct+=1
        c=0
        for reg in range(regct):
            y=align_data[reg]
            data_all[c].extend(y)
            c+=1


    BASES = ["A", "C", "G", "T"]
    BASE2INT = dict([x[::-1] for x in enumerate(BASES)])
    with open(args.report_variants,'r') as f:
        var_data: dict = {}
        header = f.readline().strip("\n").split("\t")
        header_all.extend(['strand','read_index','read_pos','template_pos', 'expected_read_base','expected_template_base','variant_read_base','variant_template_base','A1','C1','G1','T1','A2','C2','G2','T2','VarAF'])
        regct=0
        for line in f:
            x = line.strip("\n").split("\t")

            locus = x[1]
            ref_index = x[2]
            strand = x[3]
            read = x[4]
            poss = x[5:7]
            bases = x[7:11]
            onecounts = list(map(int,x[11:15]))
            twocounts = list(map(int,x[15:19]))
            if (sum(twocounts)>0) and (bases[2] in BASES):
                altbase = BASE2INT[bases[2]]
                varaf = float(twocounts[altbase])/sum(twocounts)
            else:
                varaf = 0

            if locus in var_data:
                if (var_data[locus][0] != ref_index) and (var_data[locus][2] == read): # combine counts
                    prevcounts1 = var_data[locus][9:13]
                    prevcounts2 = var_data[locus][13:17]
                    var_data[locus] = [ref_index,strand,read]
                    var_data[locus].extend(poss)
                    var_data[locus].extend(bases)
                    newcounts1 = [a+b for a,b in zip(onecounts,prevcounts1)]
                    newcounts2 = [a+b for a,b in zip(twocounts,prevcounts2)]
                    if (sum(newcounts2)>0) and (bases[2] in BASES):
                        altbase = BASE2INT[bases[2]]
                        newvaraf = float(newcounts2[altbase])/sum(newcounts2)
                    else:
                        newvaraf = 0
                    var_data[locus].extend(newcounts1)
                    var_data[locus].extend(newcounts2)
                    var_data[locus].append(newvaraf)
                # skip other read entry
            else: # first entry
                var_data[locus] = [ref_index,strand,read]
                var_data[locus].extend(poss)
                var_data[locus].extend(bases)
                var_data[locus].extend(onecounts)
                var_data[locus].extend(twocounts)
                var_data[locus].append(varaf)
                regct+=1

        c=0
        for reg in range(regct):
            if str(reg) in var_data:
                # if reads were too short to cover variant this will fail,
                # so check first
                y=var_data[str(reg)][1:]
                data_all[c].extend(y)
            else: 
                logger.info(
                    "Target variant missing from per base files %s", reg) 
            c+=1


    ########################################################################
    # Write to summary file
    output_file = args.combined_report
    with open(output_file,"w") as outfile:

        outfile.write(tabprint(header_all)+"\n")
        for data in data_all:
            outfile.write(tabprint(data)+"\n")
