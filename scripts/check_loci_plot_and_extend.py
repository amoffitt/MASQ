'''
Check in WGS BAM files for the target regions to identify any nearby SNVs or Indels
Plots the coverage and base breakdown of target regions from WGS BAMs
'''

import os
from collections import defaultdict
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np

from masq_helper_functions import setup_logger
from masq.utils.reference_genome import ReferenceGenome
from masq.utils.io import process_target_info, extend_snv_info_with_target_info
from masq.utils.io import load_snv_table, write_snv_table, tabprint
from masq.utils.seqs import BASE2INT, BASES
from masq.plots.loci_plot import loci_plot


########################################################################
# Start timer
t0 = time.time()
# Setup log file
log = setup_logger(snakemake.log,'check_loci')
log.info('Starting process')

########################################################################
# INPUT FILES AND PARAMETERS
log.info('Getting input files and parameters from snakemake object')
# WGS path to bam
# Now allows multiple BAMs
if isinstance(snakemake.input.bam,list) and isinstance(snakemake.params.wgs_name,list):
    WGS_BAM = snakemake.input.bam
    WGS_NAME = snakemake.params.wgs_name
elif isinstance(snakemake.input.bam,str) and isinstance(snakemake.params.wgs_name,str):
    WGS_BAM = [snakemake.input.bam]
    WGS_NAME = [snakemake.params.wgs_name]
else:
    log.error("WGS_BAM and WGS_BAM should be lists or strings")

sample_bams = dict(zip(WGS_NAME, WGS_BAM))

# # Input SNV table
SNV_table = open(snakemake.input.SNV_table,'r')

########################################################################
# OUTPUT FILES
log.info('Getting output files from snakemake object')
# Output plots
list_of_plot_files = snakemake.output.plots
# Updated SNP file
updated_SNV_table = snakemake.output.new_SNV_table

# Make plot folder if it doesn't exist
plotfolder = os.path.dirname(list_of_plot_files[0])
os.makedirs(plotfolder, exist_ok=True)


########################################################################
## load the genome sequence from a pickle file
log.info('Loading reference genome pickle')
# seq_pickle = snakemake.params.wgs_ref
# seq_dic = pickle.load(open(seq_pickle, 'rb'))

ref_genome = ReferenceGenome(snakemake.params.wgs_ref_fa)
ref_genome.open()
log.info('Done loading reference genome pickle')

########################################################################
# Load the input SNV table
log.info('Loading SNV table')
snv_info = load_snv_table(SNV_table)
for key,value in snv_info.items():
    log.debug(key)
    log.debug(tabprint(value))
log.info('Done loading SNV table')

########################################################################

# Process input file to extend it with strand info and get coordinates of seq
target_info = process_target_info(snv_info, ref_genome)

########################################################################
# Now go through WGS data and add non-ref bases to the input file...
log.info('Finding additional variant positions from WGS BAM')


all_targets = []
loc_ind = 0  # counter for which loci we're on

sample_bams = dict(zip(WGS_NAME, WGS_BAM))

#
# load each region from the target info read in and processed
for chrom, strand, true_start, true_end, true_targets in target_info:

    if chrom=="0": # not human / mouse sequence (GFP seq for example)
        all_targets.append([])
        image_filename = list_of_plot_files[loc_ind]
        fig = plt.figure(figsize=(20, 18))
        plt.savefig(image_filename, dpi=200, facecolor='w', edgecolor='w',
                    papertype=None, format=None,
                    transparent=False)
        plt.close()
        loc_ind += 1
    else:

        log.debug('Locus %d', (loc_ind+1))
        image_filename = list_of_plot_files[loc_ind]
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
log.info('Writing updated SNV info table')
snv_info = extend_snv_info_with_target_info(
    snv_info,
    target_info,
    all_targets,
)

with open(updated_SNV_table, 'w') as outfile:
    write_snv_table(snv_info, outfile)


########################################################################
# End timer
t1 = time.time()
td = (t1 - t0) / 60
log.info("Done in %0.2f minutes", td)

########################################################################
