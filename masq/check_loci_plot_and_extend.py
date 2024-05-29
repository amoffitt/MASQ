'''
Check in WGS BAM files for the target regions to identify any nearby SNVs or Indels
Plots the coverage and base breakdown of target regions from WGS BAMs
'''
import argparse
import logging
import sys
from typing import Optional

import matplotlib
import matplotlib.pyplot as plt

from masq.utils.reference_genome import ReferenceGenome
from masq.utils.io import process_target_info, extend_snv_info_with_target_info
from masq.utils.io import load_snv_table, write_snv_table, tabprint
from masq.plots.loci_plot import loci_plot

matplotlib.use('Agg')


logger = logging.getLogger("masq_check_loci_plot_and_extend")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="plot the coverage and base breakdown of target "
        "regions from WGS BAMs"
    )
    parser.add_argument(
        "--snv-table",
        help="SNV table filename in tab-separated format",
    )
    parser.add_argument(
        "--ref-genome",
       help="Reference genome filename in FASTA format",
    )
    parser.add_argument(
        "--wgs-bam-names",
        help="Comma separated list of WGS BAM names",
    )
    parser.add_argument(
        "--wgs-bam-files",
        help="Comma separated list of WGS BAM filenames",
    )
    parser.add_argument(
        "--output-snv-table",
        help="Output SNV table filename",
    )
    parser.add_argument(
        "plot_files",
        nargs="+",
        help="List of plot filenames",
    )
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    if argv is None:
        argv = sys.argv[1:]

    args = parse_args(argv)

    bam_names = args.wgs_bam_names.split(' ')
    bam_files = args.wgs_bam_files.split(' ')
    sample_bams = dict(zip(bam_names, bam_files))
    plot_files = args.plot_files

    with ReferenceGenome(args.ref_genome) as ref_genome:
        logger.debug("Loading SNV table")
        with open(args.snv_table, 'r') as infile:
            snv_info = load_snv_table(infile)
            for key,value in snv_info.items():
                logger.debug(key)
                logger.debug(tabprint(value))
        logger.info('Done loading SNV table')
        target_info = process_target_info(snv_info, ref_genome)


        ########################################################################
        # Now go through WGS data and add non-ref bases to the input file...
        logger.info('Finding additional variant positions from WGS BAM')

        all_targets: list[list] = []
        loc_ind = 0  # counter for which loci we're on

        #
        # load each region from the target info read in and processed
        for chrom, strand, true_start, true_end, true_targets in target_info:

            if chrom=="0": # not human / mouse sequence (GFP seq for example)
                all_targets.append([])
                image_filename = plot_files[loc_ind]
                plt.figure(figsize=(20, 18))
                plt.savefig(image_filename, dpi=200, facecolor='w', edgecolor='w',
                            papertype=None, format=None,
                            transparent=False)
                plt.close()
                loc_ind += 1
            else:

                logger.debug('Locus %d', (loc_ind+1))
                image_filename = plot_files[loc_ind]
                loc_ind += 1

                new_targets = loci_plot(
                    loc_ind,
                    image_filename,
                    chrom,
                    true_start,
                    true_end,
                    true_targets,
                    strand,
                    ref_genome,
                    sample_bams,
                )
                all_targets.append(new_targets)


        # Add new het/homo non-ref sites to original input file
        logger.info('Writing updated SNV info table')
        snv_info = extend_snv_info_with_target_info(
            snv_info,
            target_info,
            all_targets,
        )

        with open(args.output_snv_table, 'w') as outfile:
            write_snv_table(snv_info, outfile)
