'''
Checks all reads from FASTQ files for proper sequence structure
Sorts all reads by locus
Groups by varietal tag and read sequence
Allows all further processing to be done parallelized by locus
Reports on number of reads assigned to each locus
'''

import argparse
import logging
import pickle
import sys
from collections import Counter, defaultdict
from typing import Optional

import yaml
import numpy as np
import editdistance

import matplotlib.pyplot as plt

from masq.utils.io import tabprint
from masq.utils.io import load_snv_table
from masq.tags.cluster_rollup import double_counter

from masq.utils.paired_reads import PairedReads
from masq.utils.seqs import (
    hamming,
    check_tag_structure,
    convert_quality_score,
)
logger = logging.getLogger("masq_sort_data_by_tag_and_locus")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combines reports from different steps of the pipeline "
        "into one final report per sample"
    )

    parser.add_argument(
        "--fastq1",
        help="Paired-end FASTQ file 1",
    )
    parser.add_argument(
        "--fastq2",
        help="Paired-end FASTQ file 2",
    )
    parser.add_argument(
        "--snv-table",
        help="SNV table filename",
    )

    parser.add_argument(
        "--vt-counters",
        help="List of output vt_counters filenames",
    )
    parser.add_argument(
        "--vt-seq-counters",
        help="List of output vt_seq_counters filenames",
    )
    parser.add_argument(
        "--flip-counters",
        help="List of output flip_counters filenames",
    )
    parser.add_argument(
        "--ss-hamming-plot",
        help="Output hamming plot filename",
    )
    parser.add_argument(
        "--counter-report",
        help="Output counter report filename",
    )
    parser.add_argument(
        "--up2-unmatched-report",
        help="Output UP2 unmatched report filename",
    )
    parser.add_argument(
        "--ss1ss2-unmatched-report",
        help="Output SS1 and SS2 unmatched report filename",
    )
    parser.add_argument(
        "--goodtag-report",
        help="Output good tag report filename",
    )
    parser.add_argument(
        "--badtag-report",
        help="Output bad tag report filename",
    )

    parser.add_argument(
        "--configfile", default="config.yaml",
        help="Pipeline configuration YAML file",
    )
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    """Checks all reads from FASTQ files for proper sequence structure

    Sorts all reads by locus
    Groups by varietal tag and read sequence
    Allows all further processing to be done parallelized by locus
    Reports on number of reads assigned to each locus
    """
    if argv is None:
        argv = sys.argv[1:]

    args = parse_args(argv)


    ########################################################################
    # INPUT FILES AND PARAMETERS
    logger.info('Getting input files and parameters from snakemake object')
    # Input FASTQs
    fastq1 = args.fastq1
    fastq2 = args.fastq2

    with open(args.configfile, 'r') as configfile:
        config = yaml.safe_load(configfile)

        # Sequence and hamming parameters
        MAX_HAMMING_UP2 = config["UP2_hamming"]
        MAX_HAMMING_SS = config["SS_sum_hamming"]
        MIN_LEN = config["min_len"]
        TRIM_LEN = config["trim_len"]
        TAG_TEMPLATE = config["tag"]
        UP2 = config["UP2"]
        PROTOCOL = config["protocol"]
        USE_EDIT_DISTANCE = config["use_edit_distance"]

        ## compute some lengths
        UPL   = len(UP2)
        TAGL  = len(TAG_TEMPLATE)
        UPVTL = UPL + TAGL

        # Masking low quality bases
        MASK_LOWQUAL = config["mask_lowqual_bases"]
        QUAL_CUTOFF = config["qual_cutoff"]
        MAX_N_RATIO = config["max_N_ratio"]

        # Quick run parameters
        QUICK_RUN = config["quick_run"]
        QUICK_RUN_READS = config["quick_run_reads"]

    ########################################################################
    # OUTPUT FILES
    logger.info('Getting output files from snakemake object')

    # Main report
    outfile = open(args.counter_report, 'w')

    # Troubleshooting reports
    # Unmatched UP2
    up2file = open(args.up2_unmatched_report, 'w')
    # Unmatched SS1 and SS2
    ss1ss2file = open(args.ss1ss2_unmatched_report, 'w')
    # Good and bad tag structure
    goodtagfile = open(args.goodtag_report, 'w')
    badtagfile = open(args.badtag_report, 'w')

    # Plots
    # Hamming distance of SS1 and SS2
    hamming_SS_figure = args.ss_hamming_plot

    # Seq data pickles
    # Lists of files (one for each region)
    vt_counters_filelist = args.vt_counters.split(" ")
    vt_seq_counters_filelist = args.vt_seq_counters.split(" ")
    flip_counters_filelist = args.flip_counters.split(" ")



    ########################################################################
    # Load the input SNV table
    # Input SNV table
    with open(args.snv_table, 'r') as SNV_table:

        logger.info('Loading SNV table')
        snv_info = load_snv_table(SNV_table)
        for key,value in snv_info.items():
            logger.debug(key)
            logger.debug(tabprint(value))
        logger.info('Done loading SNV table')

    ########################################################################
    # Get sequence specific primer info from the SNV_table
    logger.info('Getting sequence specific primer info from table')
    SP1 = snv_info['specific-primer-1']
    SP2 = snv_info['specific-primer-2']
    ## keep track of the lengths of the primers
    ## and the longest seen of each
    SP1_len = [len(x) for x in SP1]
    SP2_len = [len(x) for x in SP2]
    MAX_SP1 = max(SP1_len)
    MAX_SP2 = max(SP2_len)

    ########################################################################
    # Setup counters for processing reads
    logger.info('Initializing counters')
    # Single counters
    counter      = 0
    good_counter = 0
    up2_counter  = 0
    flip_counter = 0
    oklength_counter = 0

    # Counters by region
    NREG = len(SP1)
    good_counter_list = [0 for _ in range(NREG+1)] # extra one for garbage
    assigned_counter_list = [0 for _ in range(NREG+1)] # extra one for garbage

    # Counter for SS1 and SS2 hamming distances
    hamming_counter = Counter()

    # Counters with compressed sequence information
    ## keep a counters for each primer pair
    ## the counters are:
    ## how many reads with each varietal tag
    ## for each varietal tag, a Counter for the associated read pairs
    ## for each varietal tag and read-pair sequence,
    ##    how often it was not-flipped (so UP2 is on read2) and flipped.
    # Lists of counters...
    vt_counters     = [Counter() for _ in SP1]
    vt_seq_counters = [defaultdict(Counter) for _ in SP1]
    flip_counters   = [defaultdict(double_counter) for _ in SP1]

    ########################################################################
    # Build read walker from 2 FASTQ files
    logger.info('Building read walker for 2 fastqs')

    ########################################################################
    # Process reads
    # Check for:
    # 1 - OK length (then trim)
    # 2 - correct UP2 sequence (then remove it and tag)
    # 3 - correct SS1 and SS2 sequences (then group reads and remove primers)
    # for regular PCR data, skip the UP2 step, check both orientations for SS1/SS2
    logger.info('Begin processing reads')
    counter = 0

    logger.debug("PROTOCOL: %s", PROTOCOL)
    if PROTOCOL != "standard PCR":
        logger.info('Assuming UP2 and VT are present (Not Standard PCR)')

    with PairedReads(fastq1, fastq2) as rwalker:
        for read_pair in rwalker:
            ## record if we flip it around so UP2 is on read2 (or SP1 is R1, SP2 is R2)
            flipped = False
            counter += 1

            # For a test run with only a small number of reads processed
            if QUICK_RUN and counter > QUICK_RUN_READS:
                logger.warning('Stopped processing reads - reached quick run cutoff')
                break

            # Write update to log file
            if counter % 10000 == 0:
                for r in range(NREG+1):
                    logger.debug(tabprint(
                    [r,
                        assigned_counter_list[r],
                        good_counter_list[r],
                        assigned_counter_list[r] / np.float64(counter),
                        good_counter_list[r] / np.float64(counter),
                        good_counter_list[r] / np.float64(assigned_counter_list[r]) ] ))
                logger.debug(' ')
                logger.debug(tabprint([counter, oklength_counter, up2_counter, good_counter, oklength_counter / np.float64(counter), up2_counter / np.float64(counter), good_counter / np.float64(counter)]))

            # Get sequences from fastq
            entry1, entry2 = read_pair
            seq1, seq2 = entry1.sequence, entry2.sequence
            qual1, qual2 = entry1.quality, entry2.quality

            # Check length
            ## if too short, skip this read pair
            if len(seq1) < MIN_LEN or len(seq2) < MIN_LEN:
                logger.debug('Skipping read for length of read pair. Sequence1: %s', seq1)
                logger.debug('Skipping read for length of read pair. Sequence2: %s', seq2)
                continue
            oklength_counter += 1

            ## otherwise, trim to length
            # Works even if TRIM_LEN is greater than length of sequence
            seq1 = seq1[:TRIM_LEN]
            seq2 = seq2[:TRIM_LEN]
            qual1 = qual1[:TRIM_LEN]
            qual2 = qual2[:TRIM_LEN]

            # check for UP2 if MASQ (not for standard PCR)
            if PROTOCOL != "standard PCR":
                # logger.info('Assuming UP2 and VT are present (Not Standard PCR)')
                ## compute hamming distance for UP2 to each read
                # check beginning of each read for match to UP2
                d1 = hamming(seq1, UP2, use_edit_distance=USE_EDIT_DISTANCE)
                d2 = hamming(seq2, UP2, use_edit_distance=USE_EDIT_DISTANCE)

                ## if UP2 matches read1, flip it
                if d1 < d2:
                    seq1, seq2 = seq2, seq1
                    qual1, qual2 = qual2, qual1
                    flipped = True
                    flip_counter += 1

                ## pull the varietal tag and the rest of seq2 (after trimming tag and UP2)
                vt         = seq2[UPL:UPVTL]
                seq2_rest  = seq2[UPVTL:]
                qual2_rest = qual2[UPVTL:]
            else: # STANDARD PCR PROTOCOL
                # fake the up2 distances so they pass check
                d1=0
                d2=100
                # fake the trimming of seq2
                seq2_rest = seq2
                qual2_rest = qual2
                # fake the varietal tag to be unique for each read
                vt = str(counter)

            ## if the match is close enough
            # if UP2 is found, move on with this read pair
            # if protocol has no UP2, move on
            if (min(d1, d2) <= MAX_HAMMING_UP2):
                # keep track of how many correct UP2s are found
                up2_counter += 1

                ## check hamming distance of seq1 to list of SP1 and seq2 to list of SP2
                ham_array1 = [hamming(seq1     , sp, use_edit_distance=USE_EDIT_DISTANCE) for sp in SP1]
                ham_array2 = [hamming(seq2_rest, sp, use_edit_distance=USE_EDIT_DISTANCE) for sp in SP2]
                ## find the "best" index for each
                ham_ind1   = np.argmin(ham_array1)
                ham_ind2   = np.argmin(ham_array2)
                ## and store the distance
                ham_dist1  = ham_array1[ham_ind1]
                ham_dist2  = ham_array2[ham_ind2]

                if PROTOCOL!="standard PCR": # only check one orientation
                    logger.debug("Regular QSD: only check one orientation")
                    ## pick the index with the smaller hamming
                    best_ind   = ham_ind1 if ham_dist1 < ham_dist2 else ham_ind2

                    ## store the value of the match into the hamming counter array
                    ## for showing size of off targetness
                    hamming_counter[(ham_array1[best_ind], ham_array2[best_ind])] += 1

                    # to continue with this read - best matches must agree and have low distances
                    hamming_sum = ham_dist1 + ham_dist2
                    hamming_match = ham_ind1 == ham_ind2
                    hamming_min = min(ham_dist1,ham_dist2)

                # For standard PCR, also check the flipped version
                else:
                    logger.debug("Standard PCR: checking flipped orientation for match")
                    ham_array3 = [hamming(seq1     , sp, use_edit_distance=USE_EDIT_DISTANCE) for sp in SP2]
                    ham_array4 = [hamming(seq2_rest, sp, use_edit_distance=USE_EDIT_DISTANCE) for sp in SP1]
                    ## find the "best" index for each
                    ham_ind3   = np.argmin(ham_array3)
                    ham_ind4   = np.argmin(ham_array4)
                    ## and store the distance
                    ham_dist3  = ham_array3[ham_ind3]
                    ham_dist4  = ham_array4[ham_ind4]

                    ## pick the index with the smaller hamming
                    min_dist=100
                    for i,ind,dist in zip([1,2,3,4],[ham_ind1,ham_ind2,ham_ind3,ham_ind4],[ham_dist1,ham_dist2,ham_dist3,ham_dist4]):
                        if dist<min_dist:
                            min_dist=dist
                            best_ind = ind
                            best_orientation = i
                    if best_orientation<3: # original orientation is good
                        logger.debug("Original Orientation")    
                        hamming_counter[(ham_array1[best_ind], ham_array2[best_ind])] += 1
                        hamming_sum = ham_dist1 + ham_dist2
                        hamming_match = ham_ind1 == ham_ind2
                        hamming_min = min(ham_dist1,ham_dist2)
                    else: # flip the reads and use the second 2 hamming comparisons
                        logger.debug("Flipped Orientation")
                        hamming_counter[(ham_array3[best_ind], ham_array4[best_ind])] += 1
                        hamming_sum = ham_dist3 + ham_dist4
                        hamming_match = ham_ind3 == ham_ind4
                        hamming_min = min(ham_dist3,ham_dist4)

                        seq1, seq2_rest = seq2_rest, seq1
                        qual1, qual2_rest = qual2_rest, qual1
                        flipped = True
                        flip_counter += 1
                        logger.debug("Flip counter: %d", flip_counter)


                ## if the sum of the distances is less than MAX Hamming and the two reads agree on the best match
                # If SP1 and SP2 look ok, move on with this read pair
                if hamming_sum <= MAX_HAMMING_SS and hamming_match:
                    ## get the rest of the reads and store the info into the counters
                    # did this read pair pass all the filters
                    good_counter += 1
                    good_counter_list[best_ind] +=1
                    # remove SP1 and SP2, take rest of read
                    ror1 = seq1     [SP1_len[best_ind]:]
                    ror2 = seq2_rest[SP2_len[best_ind]:]

                    # If mask low qual parameter setting is true, replace low quality bases with N's
                    if MASK_LOWQUAL:
                        qual1_target = qual1     [SP1_len[best_ind]:]
                        qual2_target = qual2_rest[SP2_len[best_ind]:]
                        qual1_converted = convert_quality_score(qual1_target)
                        qual2_converted = convert_quality_score(qual2_target)
                        ror1 = ''.join([ror1[x] if (qual1_converted[x]>QUAL_CUTOFF) else ('N') for x in range(len(ror1))])
                        ror2 = ''.join([ror2[x] if (qual2_converted[x]>QUAL_CUTOFF) else ('N') for x in range(len(ror2))])

                    # counter for each varietal tag, per locus
                    vt_counters[best_ind][vt] += 1
                    # store the sequences for each varietal tag, per locus
                    vt_seq_counters[best_ind][vt][(ror1, ror2)] += 1
                    # store flip status for each vt/r1/r2, per locus
                    flip_counters[best_ind][(vt, ror1, ror2)][int(flipped)] += 1
                    # this is a good read pair -  record if its tag structure is correct in "good tag" file
                    goodtagfile.write(tabprint([vt,check_tag_structure(vt,TAG_TEMPLATE)])+"\n")
                else:
                    if hamming_min<4: # ss1 and ss2 did not pass but came close
                        # Save the info for these failed reads to file
                        read_ss1= seq1[:len(SP1[best_ind])]
                        read_ss2= seq2_rest[:len(SP2[best_ind])]
                        read_ss1_ham = hamming(read_ss1,SP1[best_ind])
                        read_ss2_ham = hamming(read_ss2,SP2[best_ind])
                        read_ss1_edit = editdistance.eval(read_ss1,SP1[best_ind])
                        read_ss2_edit = editdistance.eval(read_ss2,SP2[best_ind])
                        ror2 = seq2_rest[SP2_len[best_ind]:]
                        ror1 = seq1[SP1_len[best_ind]:]
                        ss1ss2file.write(tabprint([best_ind,read_ss1,read_ss1_ham,read_ss1_edit,read_ss2,read_ss2_ham,read_ss2_edit,ror1,ror2,check_tag_structure(vt,TAG_TEMPLATE)])+"\n")
                        # this is a bad read pair - record if its tag structure is correct in "bad tag" file
                        badtagfile.write(tabprint([vt,check_tag_structure(vt,TAG_TEMPLATE)])+"\n")
                    else: # not even close to any SS's by hamming distance...
                        ss1ss2file.write(tabprint([None,seq1[:20],None,None,seq2_rest[:20],None,None,seq1[20:],seq2_rest[20:],check_tag_structure(vt,TAG_TEMPLATE)])+"\n")
                        # this is a bad read pair - record if its tag structure is correct in "bad tag" file
                        badtagfile.write(tabprint([vt,check_tag_structure(vt,TAG_TEMPLATE)])+"\n")

                # for reporting how many are assigned to each value
                if hamming_min<4: # hard coded cutoff for garbage assignments
                    assigned_counter_list[best_ind] +=1
                else:
                    assigned_counter_list[-1] += 1 # garbage bin
            else:  # UP2 did not match well on either side
                s1, s2 = entry1.sequence, entry2.sequence
                up2_seq1 = s1[:UPL]
                up2_seq2 = s2[:UPL]
                dedit1 = editdistance.eval(up2_seq1,UP2)
                dedit2 = editdistance.eval(up2_seq2,UP2)
                up2file.write(tabprint([up2_seq1,d1,dedit1,up2_seq2,d2,dedit2])+"\n")

    logger.info('Done processing reads')

    ########################################################################

    # Write to output files - counter report
    logger.info('Writing read assignment report')

    header=["Read Pairs Processed", "Pass Length Filter", "Pass Length+UP2 Filters","Pass All Filters", "Flipped Read Pair", "Fraction Length Good", "Fraction Length+UP2 Good", "Fraction All Good",  "Fraction Flipped"]
    counters=[counter, oklength_counter, up2_counter, good_counter, flip_counter, oklength_counter / np.float64(counter), up2_counter / np.float64(counter), good_counter / np.float64(counter), flip_counter/ np.float64(counter)]
    outfile.write(tabprint(header) + "\n")
    outfile.write(tabprint(counters) + "\n\n")
    outfile.write(tabprint(["Region","Assigned Reads","Pass Reads","Percent Assigned", "Percent Pass", "Percent of Assigned that Pass"])+"\n")
    for r in range(NREG+1):
        region = r if r<NREG else "NONE"
        outfile.write(tabprint(
            [region,
            assigned_counter_list[r],
            good_counter_list[r],
            assigned_counter_list[r] / np.float64(counter),
            good_counter_list[r] / np.float64(counter),
            good_counter_list[r] / np.float64(assigned_counter_list[r]) ]))
        outfile.write("\n")

    ########################################################################
    # Write condensed sequence information in counter objects to python pickle files
    # Separate data by region
    logger.info("Saving vt counters")
    for r in range(NREG):
        logger.info("Saving region %d", r)
        pickle.dump(vt_counters[r], open(vt_counters_filelist[r], 'wb'), pickle.HIGHEST_PROTOCOL)
        pickle.dump(vt_seq_counters[r], open(vt_seq_counters_filelist[r], 'wb'), pickle.HIGHEST_PROTOCOL)
        pickle.dump(flip_counters[r], open(flip_counters_filelist[r], 'wb'), pickle.HIGHEST_PROTOCOL)


    ########################################################################
    # Close output files
    logger.info('Closing output files')
    outfile.close()
    up2file.close()
    ss1ss2file.close()
    goodtagfile.close()
    badtagfile.close()

    ########################################################################
    # FIGURE - Hamming distances from SS1 and SS2
    logger.info('Plotting hamming distances from SP1 and SP2')
    N = len(SP1) # number of regions
    matrix = np.zeros(shape=(N, N), dtype=int)
    print("\n"+"Hamming distances for SS1, SS2, and counts of reads"+ "\n")
    for x in range(N):
        for y in range(N):
            print(x, y, hamming_counter[(x, y)])
            matrix[x, y] = hamming_counter[(x, y)]

    fig = plt.figure(figsize=(20, 18))
    fig.suptitle("Hamming Distance from Sequence Specific Primers", fontsize=20)
    cax = plt.imshow(np.log10(matrix), interpolation="nearest", aspect="auto")
    cbar = plt.colorbar(cax)
    cbar.set_ticks(np.arange(10))
    cbar.set_ticklabels(["$10^%d$" % d for d in np.arange(7)])
    cbar.ax.tick_params(labelsize=16)
    plt.xlabel("hamming distance from SP1 (internal)", fontsize=18)
    plt.ylabel("hamming distance from SP2 (cut-site)", fontsize=18)
    plt.tick_params(labelsize=16)
    plt.savefig(hamming_SS_figure, dpi=200, facecolor='w', edgecolor='w',
                papertype=None, format=None,
                transparent=False)
    plt.close()
