'''
Check in WGS BAM files for the target regions to identify any nearby SNVs or Indels
Plots the coverage and base breakdown of target regions from WGS BAMs
'''

import os,sys
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
import editdistance
import numpy as np
import pysam
from masq_helper_functions import tabprint
from masq_helper_functions import reverseComplement
from masq_helper_functions import convert_cigar_string
from masq_helper_functions import setup_logger, load_snv_table, write_snv_table

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
    logger.error("WGS_BAM and WGS_BAM should be lists or strings")

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
# Plotting colors
BASE_COLORS = ["#00ABBA",
               "#9ACA3C",
               "#F26421",
               "#672D8F"]

########################################################################
## load the genome sequence from a pickle file
log.info('Loading reference genome pickle')
seq_pickle = snakemake.params.wgs_ref
seq_dic = pickle.load(open(seq_pickle, 'rb'))
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
log.info('Parsing specific SNV info fields')

target_info = []
target_locs_array = [list(map(int, x.split(";") ) ) if len(x)>0 else [] for x in snv_info['target_locs']]

for i,loc in enumerate(snv_info['loc']):
    chrom=snv_info['chr'][i]
    pos = int(snv_info['posi'][i])
    target_locs = target_locs_array[i]
    aml_loc = target_locs[0]
    log.debug('Position: %d' % pos)
    log.debug('AML Loc: %d' % aml_loc)
    int_seq = snv_info['trimmed-target-seq'][i]
    length = len(int_seq)

    if chrom=="0":
        strand="+"
        start=1
        end=100
        targets=target_locs
    elif ('strand' in snv_info.keys()) and ('fragment-start' in snv_info.keys()) and ('fragment-end' in snv_info.keys()) and ('add-targets' in snv_info.keys()):
        strand=snv_info['strand'][i]
        start=snv_info['fragment-start'][i]
        end=snv_info['fragment-end'][i]
        end=snv_info['add-targets'][i]
    else:
        log.info('Inferring strand, start, and end from match to reference genome')
        # identifying start and end genomic positions of the targeted region
        # and getting target sequence from ref genome
        top_start = pos - 1 - aml_loc
        top_end = top_start + length
        top_match = seq_dic[chrom][top_start:top_end]
        # if the target sequence was on the bottom strand...
        bottom_start = pos - (length - aml_loc)
        bottom_end = bottom_start + length
        bottom_match = reverseComplement(seq_dic[chrom][bottom_start:bottom_end])

        log.debug('Target seq : %s' % int_seq)
        log.debug('Top match  : %s' % top_match)
        log.debug('Btm match  : %s' % bottom_match)

        # check which strand the sequence came from and set coordinates
        if top_match == int_seq:
            log.debug('Locus %s: positive strand' % (loc))
            strand = "+"
            start = top_start
            end = top_end
            targets = target_locs
        else:
            log.debug('Locus %s: negative strand' % (loc))
            strand = "-"
            start = bottom_start
            end = bottom_end
            targets = [length - x - 1 for x in target_locs]

    target_info.append([chrom, strand, start, end, targets])

########################################################################
# Now go through WGS data and add non-ref bases to the input file...
log.info('Finding additional variant positions from WGS BAM')
BASES = ["A", "C", "G", "T"]
#sample_names = [WGS_NAME] # EDIT TO ALLOW MULTIPLE BAMS...
BASE2INT = dict([x[::-1] for x in enumerate(BASES)])
MIN_QUAL = 20
MIN_MAP = 40
buff = 30
region_buffer = 30


all_targets = []
loc_ind = 0  # counter for which loci we're on
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

        log.debug('Locus %d' % (loc_ind+1))
        image_filename = list_of_plot_files[loc_ind]
        loc_ind += 1

        # add buffer to target sequence on either side and update values
        start = true_start - region_buffer
        end = true_end + region_buffer
        targets = [x + region_buffer for x in true_targets]
        # get local sequence
        local_seq = seq_dic[chrom][start:end]
        log.debug('Chrom: %s' % chrom)
        log.debug('start: %d' % start)
        log.debug('end: %d' % end)
        log.debug('Local_seq: %s' % local_seq)
        L = len(local_seq)
        # also in integer form
        local_int = np.array([BASE2INT[x] for x in local_seq])
        # intialize for keeping track of variants at this position
        has_something = np.zeros(L, dtype=bool)

        # for each sample
        fig = plt.figure(figsize=(20, 18))
        locus_string = "%s:%d-%d %s" % (chrom, start, end, strand)
        fig.suptitle(locus_string, fontsize=20)
        ind = 1
        for sample_name,sample_bam_fn in zip(WGS_NAME,WGS_BAM):
            log.info('Plotting locus %d' % loc_ind)
            # Set up the plot
            ax1 = fig.add_subplot(len(WGS_NAME), 1, ind)
            ax2 = ax1.twinx()
            ax1.tick_params(axis='both', labelsize=15)
            ax2.tick_params(axis='both', labelsize=15)
            ind += 1

            # open the bam file
            sample_bam = pysam.AlignmentFile(sample_bam_fn, "rb")

            # for each position of interest
            # get an iterator that walks over the reads that include our event
            read_iterator = sample_bam.fetch(
                chrom, start - buff, end + buff)  # ref has chr1,chr2,etc.

            # store as pairs
            read_names = []
            read_pair_dict = dict()
            for read in read_iterator:
                try:
                    read_pair_dict[read.qname].append(read)  # if pair is there
                except:
                    read_pair_dict[read.qname] = [read]  # if pair is absent
                    read_names.append(read.qname)  # keys to dictionary

            N = len(read_names)  # number of read pairs in this region
            # store an integer at every position for each read(pair) in interval
            # matrix (number of reads (N) by number of positions (end-start))
            read_stack = np.zeros(shape=(N, end - start), dtype=int)
            # get all read pairs
            for rcounter, read_name in enumerate(read_names):
                reads = read_pair_dict[read_name]
                for read in reads:
                    read_sequence = read.query_sequence
                    read_qualities = read.query_qualities
                    # get_aligned_pairs:
                    # returns paired list of ref seq position and read seq position
                    for rindex, pindex in read.get_aligned_pairs():
                        # rindex and pindex return ALL positions, including soft-clipped
                        # rindex: position in read (1-150 for example)
                        # pindex: position in reference genome
                        # if there's a base in the read
                        # and we are in range
                        if rindex is not None and pindex is not None:
                            # Separated these lines, added pindex None check
                            if pindex >= start and pindex < end:
                                # compute the base score
                                base_score = read_qualities[rindex]
                                # if it's good enough
                                if base_score >= MIN_QUAL:
                                    # check if there is already a value in place
                                    current = read_stack[rcounter, pindex - start]
                                    base = read_sequence[rindex]
                                    if base == "N":
                                        baseint = 0
                                    else:
                                        # A-1, C-2, G-3, T-4
                                        baseint = BASE2INT[base] + 1
                                    # if there is no value, store the value
                                    if current == 0:
                                        read_stack[rcounter, pindex - start] = baseint
                                    else:
                                        # if there is a mismatch between the two reads,
                                        # set value back to 0
                                        # this value is just for the 2 paired reads
                                        if current != baseint:
                                            read_stack[rcounter, pindex - start] = 0
            summary = []
            # iterating over numpy array - by rows
            # transpose first to iterate over positions
            for x in read_stack.transpose():
                # gets counts of N,A,C,G,T per position as array
                # append to summary list
                summary.append(np.bincount(x, minlength=5))
            # convert summary to array, fills by row
            # drop the N count, and transpose again
            # .T tranposes, only if num dimensions>1
            summary = np.array(summary)[:, 1:].T
            # now we have base (4) by position array as summary
            # base_cover: coverage at each position (sum of A/C/G/T counts)
            base_cover = np.sum(summary, axis=0)
            # base_ratio: A/C/G/T counts over total coverage
            # aka frequency of each base
            base_ratio = summary.astype(float) / np.maximum(1, base_cover)
            # update has something
            # EDIT - 0 coverage is not has_something??
            has_something += ( (base_ratio[local_int, np.arange(L)] < 0.9) & (base_cover>3) ) # reference ratio is less than 0.9

            ########################################################################
            # Plot variants in each region

            # plot the coverage first
            ax1.plot(base_cover, color='k', label='coverage', alpha=0.5)
            # draw lines for boundaries of event
            # and targets
            if strand == "+":
                ax2.axvline(region_buffer-0.5, color='g', lw=2)
                ax2.axvline(L-region_buffer-0.5, color='g', linestyle='--', lw=2)
            else:
                ax2.axvline(region_buffer-0.5, color='g', linestyle='--', lw=2)
                ax2.axvline(L-region_buffer-0.5, color='g', lw=2)
            for pos in targets:
                ax1.axvline(pos, color='r', lw=2)

            # Plot colored circle for each base at position vs. frequency
            for BASE_COLOR, BASE, yvals in zip(BASE_COLORS, BASES, base_ratio):
                ax2.plot(yvals, 'o', markersize=13, label=BASE, color=BASE_COLOR)
                base_filter = local_int == BASE2INT[BASE]  # same as ref base
                x = np.arange(L)
                # label ref bases with  black circle
                ax2.plot(x[base_filter], yvals[base_filter], 'o',
                         markersize=13, mfc="None", mec='black', mew=2)
            ax1.set_ylabel(sample_name, fontsize=25)

            # for non-ref sites...
            # label number of reads for non-ref base in red
            for hs_ind in np.where(has_something)[0]:
                base_counts, total = summary[:, hs_ind], base_cover[hs_ind]
                ref_base = local_int[hs_ind]
                for bind in range(4):
                    # loop through A/C/G/T, find the non-zero counts
                    if bind != ref_base and base_counts[bind] > 0:
                        # text_string = r"$\frac{%d}{%d}$" % (
                            # base_counts[bind],total)
                        # text_string = r"%d/%d" % (base_counts[bind], total)
                        text_string = r"%d" % (base_counts[bind])
                        # add the count to the plot
                        ax2.text(hs_ind + 1, base_ratio[bind, hs_ind], text_string,
                                 fontsize=20, color='red', weight='semibold')
            ax2.set_xlim(-10, 10+L)
            ax2.set_ylim(0, 1)
            ax2.legend(loc="upper center", numpoints=1, ncol=4,
                       fontsize=18, columnspacing=0.8, handletextpad=-0.2)
            ax1.legend(loc="lower center", numpoints=2,
                       fontsize=18, columnspacing=0.5)

        plt.tight_layout(rect=(0, 0, 1, 0.98))
        plt.savefig(image_filename, dpi=200, facecolor='w', edgecolor='w',
                    papertype=None, format=None,
                    transparent=False)
        plt.close()

    ########################################################################
        # Add non-ref het/homo sites to "target" list
        rlen = true_end - true_start
        new_targets = np.where(has_something)[0] - region_buffer
        if strand == "-":
            new_targets = true_end - true_start - new_targets - 1

        # new_targets>=0 (in original target region)
        new_targets = new_targets[(new_targets >= 0) * (new_targets < rlen)]
        all_targets.append(new_targets)

    ########################################################################

# Add new het/homo non-ref sites to original input file
log.info('Writing updated SNV info table')
outfile = open(updated_SNV_table, 'w')

snv_info['add-targets']=[]
snv_info['strand']=[]
snv_info['fragment-start']=[]
snv_info['fragment-end']=[]
print("Start loop")
for i, (new_targets, more_info) in enumerate(zip(all_targets, target_info)):
    print(i)
    prev_targets = list(map(int, snv_info['target_locs'][i].split(";")))
    add_targets = [x for x in new_targets if x not in prev_targets]
    snv_info['add-targets'].append(";".join(list(map(str, add_targets))))
    snv_info['strand'].append(more_info[1])
    snv_info['fragment-start'].append(more_info[2])
    snv_info['fragment-end'].append(more_info[3])
print("End loop")
print(snv_info)

write_snv_table(snv_info,outfile)

########################################################################
# Close files
outfile.close()

########################################################################
# End timer
t1 = time.time()
td = (t1 - t0) / 60
log.info("Done in %0.2f minutes" % td)

########################################################################
