'''
Consensus base calling and reporting of template counts
Reports give template counts for different filters
Exact and greater than X reads per template
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
    reverse_complement,
    test_pair,
    test_pair_withNs,
)
logger = logging.getLogger("masq_all_base_report")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Consensus base calling and reporting of template counts. "
        "Reports give template counts for different filters. "
        "Exact and greater than X reads per template"
    )

    parser.add_argument(
        "--vt-counter",
        help="vt_counter filename",
    )
    parser.add_argument(
        "--vt-seq-counter",
        help="vt_seq_counter filename",
    )
    parser.add_argument(
        "--flip-counter",
        help="flip_counter filename",
    )
    parser.add_argument(
        "--base-count-report",
        help="Output base count report filename",
    )
    parser.add_argument(
        "--base-count-report-per-base",
        help="Output base count report per base filename",
    )
    parser.add_argument(
        "--alignment-report",
        help="Output alignment report filename",
    )
    parser.add_argument(
        "--within-tag-errors",
        help="Output within tag errors filename",
    )
    parser.add_argument(
        "--within-tag-errors-table",
        help="Output within tag errors table filename",
    )
    parser.add_argument(
        "--unaligned-reads",
        help="Output unaligned reads filename",
    )
    parser.add_argument(
        "--snv-table",
        help="Input SNV table filename",
    )
    parser.add_argument(
        "--region",
        help="Region index",
    )

    parser.add_argument(
        "--configfile", default="config.yaml",
        help="Pipeline configuration YAML file",
    )
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    """Consensus base calling and reporting of template counts.

    Reports give template counts for different filters
    Exact and greater than X reads per template
    """
    if argv is None:
        argv = sys.argv[1:]

    args = parse_args(argv)

    ########################################################################
    # INPUT FILES AND PARAMETERS
    logger.info('Getting input files and parameters from snakemake object')

    # Read data stored in counters
    vt_counter_filename = args.vt_counter
    vt_seq_counter_filename = args.vt_seq_counter

    # Input SNV table
    SNV_table = open(args.snv_table,'r')

    # WHICH REGION ARE WE PROCESSING
    REGION = args.region

    with open(args.configfile, 'r') as configfile:
        config = yaml.safe_load(configfile)


        # Sequence and hamming parameters
        MAX_HAMMING_TARGET = config["target_hamming"]
        # Possible minimum coverage cutoffs
        # need to see this many reads at that position to call a base for each tag
        MIN_COVERS = config["coverage_list"]
        # Error rate allowed within reads of the same tag
        MAX_ERROR = config["base_error_rate"]
        CUTOFF = (1. - MAX_ERROR)

        # Masking low quality bases
        MASK_LOWQUAL = config["mask_lowqual_bases"]
        MAX_N_RATIO = config["max_N_ratio"]
        FILTER_NS = config["filter_ns"]
        trim_len = config["trim_len"]
        filter_ns = config["filter_ns"]
        UP2 = config["UP2"]
        tag = config["tag"]


    ########################################################################
    EXACT_COVERS = np.arange(1,51)

    ########################################################################
    POSSIBLE_PROPORTIONS = []
    PROPORTIONS_BY_COV = dict()
    for coverage in MIN_COVERS:
        a=np.linspace(0,1,coverage+1)
        POSSIBLE_PROPORTIONS.extend(a)
        PROPORTIONS_BY_COV[coverage]=a
    POSSIBLE_PROPORTIONS = np.unique(POSSIBLE_PROPORTIONS)

    ########################################################################
    # OUTPUT FILES
    logger.info('Getting output files from snakemake object')

    # All base report for specific region
    # Original report style
    base_count_per_pos_file = open(args.base_count_report, 'w')
    # One base per line (with >=X reads per tag, and =X reads per tag)
    base_count_per_base_file = open(args.base_count_report_per_base, 'w')

    # Alignment counters report for region
    alignment_rate_file = open(args.alignment_report, 'w')

    # Unaligned reads
    unaligned_reads_file = open(args.unaligned_reads, 'w')

    # Within tag error
    withintagerrors_file = open(args.within_tag_errors_table, 'w')


    ########################################################################
    # Nucleotide <-> Integer mapping
    NUCS = ["A", "C", "G", "T"]
    BASES = ["A", "C", "G", "T", "N"]
    BASE2INT = dict(x[::-1] for x in enumerate(BASES))

    # FOR WITHIN TAG ERROR
    REF_TRINUCS = [x+y+z for x in NUCS for y in NUCS for z in NUCS]
    # All pairs of REFS and ALT trinucs
    # PAIRED_TRINUCS=[ (a , a[0]+x+a[2], b) for b in ('R1','R2') for a in REF_TRINUCS for x in NUCS if x!=a[1] ]
    PAIRED_TRINUCS=[ (a , a[0]+x+a[2]) for a in REF_TRINUCS for x in NUCS if x!=a[1] ]
    # Trinuc to Index
    TRINUC2INT = dict(x[::-1] for x in enumerate(PAIRED_TRINUCS)) # 192
    # Array of error counts
    ERR_RANGE1 = np.arange(-0.1,1.0,0.1)
    ERR_RANGE2 = np.arange(0,1.01,0.1)
    s1 = len(TRINUC2INT)
    s2 = len(ERR_RANGE1)
    s3 = len(('R1','R2'))
    READS_PER_TAG_CUTOFF = 10

    ########################################################################
    # Load sequence data
    logger.info('Loading read data from pickles')
    vt_counter = pickle.load(open(vt_counter_filename, 'rb')) # only one region
    vt_seq_counter = pickle.load(open(vt_seq_counter_filename, 'rb'))
    logger.info('Done loading read data from pickles')

    ########################################################################
    # Write headers to output files
    logger.info('Writing headers to file')
    # Alignment counters report for region
    header=["Region Index","Ref Sequence Index","Number of Tags", "Number of Aligned Tags", "Fraction Aligned Tags","Number of Processed Reads", "Number of Aligned Reads", "Fraction Aligned Reads"]
    alignment_rate_file.write(tabprint(header) + "\n")

    # Original report style
    position_header = ['target_base', 'locus_index', 'reference_index', 'target_strand', 'read_index', 'read_pos', 'template_pos', 'expected_read_base', 'expected_template_base', 'variant_read_base', 'variant_template_base']
    hlists = [["".join([nuc, str(cover)]) for nuc in NUCS] for cover in MIN_COVERS]
    header = position_header + [item for sublist in hlists for item in sublist]
    base_count_per_pos_file.write(tabprint(header) + "\n")

    # One base per line (with >=X reads per tag, and =X reads per tag)
    position_header = ['target_base', 'locus_index', 'reference_index', 'target_strand', 'read_index', 'read_pos', 'template_pos', 'expected_read_base', 'expected_template_base', 'variant_read_base', 'variant_template_base', 'read_ref_trinuc','template_ref_trinuc']
    base_header = ['read_alt_trinuc','template_alt_trinuc','read_alt_base','proportion_of_RPT','total_count','alt_count','at_least_x_reads','exactly_x_reads','proportions_exactly_x_reads']
    header = position_header + base_header
    base_count_per_base_file.write(tabprint(header) + "\n")


    ########################################################################
    # Load the input SNV table
    logger.info('Loading SNV table')
    snv_info = load_snv_table(SNV_table)
    for key,value in snv_info.items():
        logger.debug(key)
        logger.debug(tabprint(value))
    logger.info('Done loading SNV table')

    logger.info('Parsing specific SNV info fields')
    SP1 = snv_info['specific-primer-1']
    SP2 = snv_info['specific-primer-2']
    region_array = snv_info['trimmed-target-seq'] # target region seqeunce
    target_locs_array = [list(map(int, x.split(";") ) ) if len(x)>0 else [] for x in snv_info['target_locs']]  # Store AML loci in target_locs_array
    add_locs_array = [list(map(int, x.split(";") ) ) if len(x)>0 else [] for x in snv_info['add-targets']]  # Store other loci in add_locs_array
    region_strand_array = snv_info['strand'] # strand of target region sequence
    refbase = [x[0] for x in snv_info['ref-alt_allele']]
    altbase = [x[2] for x in snv_info['ref-alt_allele']]

    ########################################################################
    # Process region specific information, including indels
    # Store multiple reference sequences and positions resulting from indels in list
    logger.info('Processing region specific info')
    # Only processing one specific region in this script: REGION
    region_index = int(REGION)
    ## load region specific information
    region_strand = region_strand_array[region_index]
    var_base = altbase[region_index]

    # Now things that vary with indels
    list_region_seq = [region_array[region_index]]
    list_region_rc  = [reverse_complement(region_array[region_index])]
    list_region_len = [len(region_array[region_index])]

    ## get the AML target positions (may be more than one)
    # add additional loci to target list
    ## and their index for both forward (read1) and reverse (read2)
    targets    = target_locs_array[region_index] + add_locs_array[region_index]
    targets_rc = [len(list_region_seq[0]) - pos - 1 for pos in targets]
    ## aml only targets
    aml_only    = target_locs_array[region_index]
    aml_only_rc = [len(list_region_seq[0]) - pos - 1 for pos in aml_only]
    # so we can add another entry if there is an indel
    list_targets = [targets]
    list_targets_rc = [targets_rc]
    list_aml_only = [aml_only]
    list_aml_only_rc = [aml_only_rc]

    if 'indel_start' in snv_info:
        if len(snv_info['indel_start'][region_index])>0:
            logger.warning('Indels included for this region. Assuming only one indel is present')

            istart = int(snv_info['indel_start'][region_index])
            ilen = int(snv_info['indel_length'][region_index])
            iseq = snv_info['indel_seq'][region_index]

            new_region_seq = list_region_seq[0]
            if ilen<0: # deletion
                new_region_seq = new_region_seq[0:istart] + new_region_seq[(istart-ilen):]
            else: # insertion
                new_region_seq = new_region_seq[0:istart] + iseq + new_region_seq[istart:]

            list_region_seq.append(new_region_seq)
            list_region_rc.append(reverse_complement(new_region_seq))
            list_region_len.append(len(new_region_seq))

            # Adjust targets
            new_targets = [x if x<istart else (x+ilen) for x in targets]
            new_targets_rc = [len(new_region_seq) - pos - 1 for pos in new_targets]
            new_aml_only = [x if x<istart else (x+ilen) for x in aml_only]
            new_aml_only_rc = [len(new_region_seq) - pos - 1 for pos in new_aml_only]
            list_targets.append(new_targets)
            list_targets_rc.append(new_targets_rc)
            list_aml_only.append(new_aml_only)
            list_aml_only_rc.append(new_aml_only_rc)

    ########################################################################
    # Output variables
    list_WITHIN_TAG_ERRS = []

    # Loop over each possible reference sequence and go through tags to align
    logger.info('Parsing each tag and associated reads')
    for refiter in range(len(list_targets)):
        logger.info('Target sequence iteration: %d' % refiter)
        logger.info('Parsing %d tags in vt_counter' % len(vt_counter))
        logger.info('Parsing %d tags in vt_seq_counter' % len(vt_seq_counter))

        # Extract info for this reference sequence
        logger.info('Extracting info specific to this reference sequence')
        region_seq = list_region_seq[refiter]
        logger.info('Reference sequence: %s' % region_seq)
        region_rc = list_region_rc[refiter]
        region_len = list_region_len[refiter]
        targets = list_targets[refiter]
        targets_rc = list_targets_rc[refiter]
        aml_only = list_aml_only[refiter]
        aml_only_rc = list_aml_only_rc[refiter]

        # Max sequence lengths to count bases on
        # R1: min(regionlen, trimlen - SP1len)
        # R2: min(regionlen, trimlen - SP2len - UP2len - Taglen)
        L1 = min(region_len, trim_len - len(SP1[region_index]))
        L2 = min(region_len, trim_len - len(SP2[region_index]) - len(UP2) - len(tag))
        logger.debug('L1 is %d' % L1)
        logger.debug('L2 is %d' % L2)
        ## and establish default ranges
        range_L1 = np.arange(L1)
        range_L2 = np.arange(L2)

        ## and figure out which the AML target positions are active for each read
        r1_targets = [target for target in targets if target < L1]
        r2_targets = [target for target in targets_rc if target < L2]

        r1_aml_only = [target for target in aml_only if target < L1]
        r2_aml_only = [target for target in aml_only_rc if target < L2]

        # Setup counters
        logger.info('Setting up counters')
        # Basic alignment counters
        tag_counter = 0
        all_counter = 0
        hit_counter = 0
        aligned_tag_counter = 0
        # Within tag errors
        WITHIN_TAG_ERRS=np.zeros(shape=(s1,s2,s3), dtype=int)
        # store count data at each position subject to the coverage constraints in MIN_COVERS
        # read pos x coverage values x base
        base_counter1 = np.zeros(shape=(L1, len(MIN_COVERS), 4), dtype=int)
        base_counter2 = np.zeros(shape=(L2, len(MIN_COVERS), 4), dtype=int)

        base_counter1_exact = np.zeros(shape=(L1, len(EXACT_COVERS), 4), dtype=int) # changed to include all values 1 to 50
        base_counter2_exact = np.zeros(shape=(L2, len(EXACT_COVERS), 4), dtype=int)

        base_counter1_proportions = np.zeros(shape=(L1, len(MIN_COVERS), 4, len(POSSIBLE_PROPORTIONS)), dtype=int)
        base_counter2_proportions = np.zeros(shape=(L2, len(MIN_COVERS), 4, len(POSSIBLE_PROPORTIONS)), dtype=int)

        ########################################################################

        ## for each tag
        for tag, tag_count in vt_counter.most_common():
            tag_counter += 1
            tag_aligned = False
            RPT = 0
            if tag_counter % 10000 == 0:
                logger.info(tabprint([tag_counter,
                                np.max(base_counter1),
                                np.max(base_counter2),
                                all_counter,
                                hit_counter,
                                float(hit_counter) / max(1,all_counter)]))
            ## get sequence counter associated with the tag
            seq_counter = vt_seq_counter[tag]
            ## record aggregate base call data (tag-data-1 and -2)
            td1    = np.zeros(shape=(5, L1), dtype=int)
            td2    = np.zeros(shape=(5, L2), dtype=int)
            ## iterate through read pairs (r1, r2)
            for (r1, r2), seq_count in seq_counter.most_common():
                RPT += seq_count
                ## trim read pairs as needed to fit target region
                r1 = r1[:L1]
                r2 = r2[:L2]
                aligned=False
                # If reads are shorter than expected, add N's...
                if len(r1)<L1:
                    R1_long=r1+''.join(['N' for z in range(len(r1),L1)])
                    logger.debug("R1 length modified: %s" % R1_long)
                    r1=R1_long
                if len(r2)<L2:
                    R2_long=r2+''.join(['N' for z in range(len(r2),L2)])
                    logger.debug("R2 length modified: %s" % R2_long)
                    r2=R2_long

                all_counter += seq_count # number of reads seen


                ## measure the match to reference sequences
                if MASK_LOWQUAL or FILTER_NS:
                    score = test_pair_withNs(r1, r2, region_seq, region_rc, maxNratio=MAX_N_RATIO)
                    ## omit the target positions from the score
                    # N's are not penalized so don't fix them
                    for pos in r1_targets:
                        if (pos<len(r1)):
                            score -= int((r1[pos] != region_seq[pos]) and (r1[pos] != 'N'))
                    for pos in r2_targets:
                        if (pos<len(r2)):
                            score -= int((r2[pos] != region_rc [pos]) and (r2[pos] != 'N'))
                else:
                    score = test_pair(r1, r2, region_seq, region_rc)
                    ## omit the target positions from the score
                    for pos in r1_targets:
                        if (pos<len(r1)):
                            score -= int(r1[pos] != region_seq[pos])
                    for pos in r2_targets:
                        if (pos<len(r2)):
                            score -= int(r2[pos] != region_rc [pos])


                if score  <= MAX_HAMMING_TARGET:
                    hit_counter += seq_count # reads match sequence of interest!
                    r1_int = [BASE2INT[x] for x in r1]
                    r2_int = [BASE2INT[x] for x in r2]
                    ## and add to the base call data according to the count
                    np.add.at(td1, [r1_int, range_L1], seq_count)
                    np.add.at(td2, [r2_int, range_L2], seq_count)
                    aligned=True
                    tag_aligned=True

                else: # If Hamming Distance wasn't a good match, try a 1-2 base shift in the sequence
                    origscore = score
                    r1_shifted = 'N' + r1[:-1] # shift one to the right
                    r2_shifted = 'N' + r2[:-1] # shift one to the right
                    score = test_pair_withNs(r1_shifted, r2_shifted, region_seq, region_rc, maxNratio=MAX_N_RATIO)
                    for pos in r1_targets:
                        if (pos<len(r1)):
                            score -= int((r1_shifted[pos] != region_seq[pos]) and (r1_shifted[pos] != 'N'))
                    for pos in r2_targets:
                        if (pos<len(r2)):
                            score -= int((r2_shifted[pos] != region_rc [pos]) and (r2_shifted[pos] != 'N'))
                    if score  <= MAX_HAMMING_TARGET: # Now, after shifting, if the score is good enough...
                        logger.debug("Left shift 1 worked")
                        hit_counter += seq_count # reads match sequence of interest!
                        r1_int = [BASE2INT[x] for x in r1_shifted]
                        r2_int = [BASE2INT[x] for x in r2_shifted]
                        ## and add to the base call data according to the count
                        np.add.at(td1, [r1_int, range_L1], seq_count) # problem here
                        np.add.at(td2, [r2_int, range_L2], seq_count)
                        aligned=True
                        tag_aligned=True
                    else:
                        r1_shifted = r1[1:] + 'N' # shift one to the left
                        r2_shifted = r2[1:] + 'N' # shift one to the left
                        score = test_pair_withNs(r1_shifted, r2_shifted, region_seq, region_rc, maxNratio=MAX_N_RATIO)
                        for pos in r1_targets:
                            if (pos<len(r1)):
                                score -= int((r1_shifted[pos] != region_seq[pos]) and (r1_shifted[pos] != 'N'))
                        for pos in r2_targets:
                            if (pos<len(r2)):
                                score -= int((r2_shifted[pos] != region_rc [pos]) and (r2_shifted[pos] != 'N'))
                        if score  <= MAX_HAMMING_TARGET: # Now, after shifting, if the score is good enough...
                            logger.debug("Right shift 1 worked")
                            hit_counter += seq_count # reads match sequence of interest!
                            r1_int = [BASE2INT[x] for x in r1_shifted]
                            r2_int = [BASE2INT[x] for x in r2_shifted]
                            ## and add to the base call data according to the count
                            np.add.at(td1, [r1_int, range_L1], seq_count)
                            np.add.at(td2, [r2_int, range_L2], seq_count)
                            aligned=True
                            tag_aligned=True

                # Write out unaligned reads
                if not aligned:
                    unaligned_reads_file.write(str(REGION)+"\t"+r1+"\t"+r2+"\t"+"\t"+str(origscore)+"\t"+str(score)+"\t"+tag+"\n")

            if tag_aligned:
                aligned_tag_counter +=1
            ## now compute some statistics over the base call information
            total1   = np.sum(td1[:4], axis = 0)    # total coverage (r1) at each position
            max_ind1 = np.argmax(td1, axis=0)       # index of maximal base (r1)
            max_val1 = np.max(td1, axis=0)          # count value at maximal base (r1)
            total2   = np.sum(td2[:4], axis = 0)    # as above for r2.
            max_ind2 = np.argmax(td2, axis=0)
            max_val2 = np.max(td2, axis=0)

            logger.debug(["tagalign",tag,str(tag_aligned),str(RPT)])

            ## for each coverage cut off
            for cind, EXACT_COVER in enumerate(EXACT_COVERS):

                # Also want to look at counts for exactly x reads per tag        
                if EXACT_COVER==EXACT_COVERS[-1]: # catch all for 50 or more
                    strong1  = (total1 >= EXACT_COVER) * (max_val1 >= CUTOFF*total1) * (max_ind1 != 4)
                    strong2  = (total2 >= EXACT_COVER) * (max_val2 >= CUTOFF*total2) * (max_ind2 != 4)
                else:
                    strong1  = (total1 == EXACT_COVER) * (max_val1 >= CUTOFF*total1) * (max_ind1 != 4)
                    strong2  = (total2 == EXACT_COVER) * (max_val2 >= CUTOFF*total2) * (max_ind2 != 4)

                base_counter1_exact[strong1, cind, max_ind1[strong1]] += 1
                base_counter2_exact[strong2, cind, max_ind2[strong2]] += 1


            for cind, MIN_COVER in enumerate(MIN_COVERS):
                ## test against coverage and base_ratio cutoffs
                # Check that max index is not an N
                strong1  = (total1 >= MIN_COVER) * (max_val1 >= CUTOFF*total1) * (max_ind1 != 4)
                strong2  = (total2 >= MIN_COVER) * (max_val2 >= CUTOFF*total2) * (max_ind2 != 4)
                ## increment those positions that pass the tests
                base_counter1[strong1, cind, max_ind1[strong1]] += 1
                base_counter2[strong2, cind, max_ind2[strong2]] += 1

                # Also want to keep track of bases that don't pass the 0.8 cutoff
                strong1  = (total1 == MIN_COVER) * (max_ind1 != 4)
                strong2  = (total2 == MIN_COVER) * (max_ind2 != 4)

                for B in range(0,4):
                    p_obs1 = np.true_divide(td1[B,:],total1)
                    p_obs2 = np.true_divide(td2[B,:],total2)

                    for pind, PROPORTION in enumerate(POSSIBLE_PROPORTIONS):
                        if PROPORTION in PROPORTIONS_BY_COV[MIN_COVER]:

                                strong1p = (strong1) * (p_obs1==PROPORTION)
                                strong2p = (strong2) * (p_obs2==PROPORTION)

                                base_counter1_proportions[strong1p, cind, B, pind] += 1
                                base_counter2_proportions[strong2p, cind, B, pind] += 1


            #######################################################################
            # WITHIN TAG ERROR
            # For R1, use td1, region_seq
            # For R2, use td2, region_rc
            # Positions to skip: r1_targets, r2_targets

            if tag_count > READS_PER_TAG_CUTOFF:
                # R1
                fractions = td1[:4,]/total1 # Get A/C/G/T fractions from counts
                errind = np.ceil(fractions*10) # gets index for the 10 err bins
                for pos in range(1,(L1-1)): # skip first and last to get trinucs
                    if pos not in r1_targets:
                        ref_base= region_seq[pos]
                        ref_trinuc= region_seq[(pos-1):(pos+2)]
                        alt_bases = [x for x in NUCS if x!=ref_base]

                        for alt in alt_bases:
                            alt_trinuc = ref_trinuc[0] + alt + ref_trinuc[2]
                            try:
                                e = int(errind[BASE2INT[alt],pos])
                            except:
                                e = None
                            if e is not None:
                                WITHIN_TAG_ERRS[TRINUC2INT[(ref_trinuc,alt_trinuc)],e,0] += 1


                # R2
                fractions = td2[:4,]/total2 # Get A/C/G/T fractions from counts
                errind = np.ceil(fractions*10) # gets index for the 10 err bins
                for pos in range(1,(L2-1)): # skip first and last to get trinucs
                    if pos not in r2_targets:
                        ref_base= region_rc[pos]
                        ref_trinuc= region_rc[(pos-1):(pos+2)]
                        alt_bases = [x for x in NUCS if x!=ref_base]

                        for alt in alt_bases:
                            alt_trinuc = ref_trinuc[0] + alt + ref_trinuc[2]
                            try:
                                e = int(errind[BASE2INT[alt],pos])
                            except:
                                e = None
                            if e is not None:
                                WITHIN_TAG_ERRS[TRINUC2INT[(ref_trinuc,alt_trinuc)],e,1] += 1


        #######################################################################
        logger.info('Done with all tags for this reference sequence')
        list_WITHIN_TAG_ERRS.append(WITHIN_TAG_ERRS)

        # Write counters results to file
        counters=[region_index, refiter, tag_counter, aligned_tag_counter, aligned_tag_counter / float(max(1,tag_counter)), all_counter, hit_counter, hit_counter / float(max(1,all_counter))]
        alignment_rate_file.write(tabprint(counters) + "\n")

        ########################################################################
        # ORIGINAL REPORT FORMAT

            ## write the output for each position over all coverage levels
        for pos, data in enumerate(base_counter1):
            # indicator for what type of variant position it is
            # 0 - ref base
            # 1 - non-aml variant base
            # 2 - AML base
            aml_loc = int(pos in r1_targets) + int(pos in r1_aml_only)
            if aml_loc==2:
                variant_read_base = var_base
                variant_template_base = var_base
            else:
                variant_read_base = '-'
                variant_template_base = '-'
            pos_in_region = pos
            outline = tabprint([aml_loc,
                                region_index,
                                refiter,
                                region_strand,
                                1,
                                pos,
                                pos_in_region,
                                region_seq[pos],
                                region_seq[pos_in_region],
                                variant_read_base,
                                variant_template_base] + list(data.ravel()))
            base_count_per_pos_file.write(outline)
            base_count_per_pos_file.write("\n")

        for pos, data in enumerate(base_counter2):
            aml_loc = int(pos in r2_targets) + int(pos in r2_aml_only)
            if aml_loc==2:
                variant_read_base = reverse_complement(var_base)
                variant_template_base = var_base
            else:
                variant_read_base = '-'
                variant_template_base = '-'
            pos_in_region = region_len - pos - 1
            outline = tabprint([aml_loc,
                                region_index,
                                refiter,
                                region_strand,
                                2,
                                pos,
                                pos_in_region,
                                region_rc[pos],
                                region_seq[pos_in_region],
                                variant_read_base,
                                variant_template_base] + list(data.ravel()))
            base_count_per_pos_file.write(outline)
            base_count_per_pos_file.write("\n")

        ########################################################################
        # NEW REPORT FORMAT
        ## write the output for each position over all coverage levels
        ## Separate one line per base and per coverage level
        ## Do for exactly and at least X reads per base

        ########################################################################
        # counts by proportion of tag
        print("BC1_proportions")
        print(base_counter1_proportions.shape)
        for pos, data in enumerate(base_counter1_proportions):
            aml_loc = int(pos in r1_targets) + int(pos in r1_aml_only)
            if aml_loc==2:
                variant_read_base = var_base
                variant_template_base = var_base
            else:
                variant_read_base = '-'
                variant_template_base = '-'
            pos_in_region = pos

            if pos==0:
                read_trinuc='X'+region_seq[0:(pos+2)]
                template_trinuc='X'+region_seq[0:(pos_in_region+2)]
            elif pos==(len(region_seq)-1):
                read_trinuc=region_seq[(pos-1):(pos+2)]+'X'
                template_trinuc=region_seq[(pos_in_region-1):(pos_in_region+2)]+'X'
            else:
                read_trinuc=region_seq[(pos-1):(pos+2)]
                template_trinuc=region_seq[(pos_in_region-1):(pos_in_region+2)]

            print("Position: %d: read_trinuc: %s" % (pos,read_trinuc))
            print("Position: %d: template_trinuc: %s" % (pos,template_trinuc))

            for cindex,cov in enumerate(MIN_COVERS):
                for nindex,nuc in enumerate(NUCS):
                    for pind, PROPORTION in enumerate(POSSIBLE_PROPORTIONS):
                        if PROPORTION in PROPORTIONS_BY_COV[cov]:

                            count = data[(cindex,nindex,pind)]
                            total = np.sum(data[cindex,0]) # sum up proportions for A (all tags are either 100% A, 50% A or 0% A for 2 RPT)

                            alt_read_trinuc = read_trinuc[0] + nuc + read_trinuc[2]
                            alt_template_trinuc = template_trinuc[0] + nuc + template_trinuc[2]

                            baseinfo = [alt_read_trinuc,alt_template_trinuc,nuc,PROPORTION,total,count,'.','.',cov] 

                            outline = tabprint([aml_loc,
                                                region_index,
                                                refiter,
                                                region_strand,
                                                1,
                                                pos,
                                                pos_in_region,
                                                region_seq[pos], # read base
                                                region_seq[pos_in_region], # template base
                                                variant_read_base,
                                                variant_template_base,
                                                read_trinuc,
                                                template_trinuc] +
                                                baseinfo)
                            base_count_per_base_file.write(outline)
                            base_count_per_base_file.write("\n")
        ##########################################
        print("BC2_proportions")
        print(base_counter2_proportions.shape)
        for pos, data in enumerate(base_counter2_proportions):
            aml_loc = int(pos in r2_targets) + int(pos in r2_aml_only)
            if aml_loc==2:
                variant_read_base = reverse_complement(var_base)
                variant_template_base = var_base
            else:
                variant_read_base = '-'
                variant_template_base = '-'
            pos_in_region = region_len - pos - 1

            if pos==0:
                read_trinuc='X'+region_rc[0:(pos+2)]
                template_trinuc=region_seq[(pos_in_region-1):(pos_in_region+2)]+'X'
            elif pos==(len(region_seq)-1):
                read_trinuc=region_rc[(pos-1):(pos+2)]+'X'
                template_trinuc='X'+region_seq[0:(pos_in_region+2)]
            else:
                read_trinuc=region_rc[(pos-1):(pos+2)]
                template_trinuc=region_seq[(pos_in_region-1):(pos_in_region+2)]
                print("Position: %d: read_trinuc: %s" % (pos,read_trinuc))
                print("Position: %d: template_trinuc: %s" % (pos,template_trinuc))

            for cindex,cov in enumerate(MIN_COVERS):
                for nindex,nuc in enumerate(NUCS):
                    for pind, PROPORTION in enumerate(POSSIBLE_PROPORTIONS):
                        if PROPORTION in PROPORTIONS_BY_COV[cov]:

                            count = data[(cindex,nindex,pind)]
                            total = np.sum(data[cindex,0])

                            alt_read_trinuc = read_trinuc[0] + nuc + read_trinuc[2]
                            alt_template_trinuc = template_trinuc[0] + reverse_complement(nuc) + template_trinuc[2]

                            baseinfo = [alt_read_trinuc,alt_template_trinuc,nuc,PROPORTION,total,count,'.','.',cov] 

                            outline = tabprint([aml_loc,
                                                region_index,
                                                refiter,
                                                region_strand,
                                                2,
                                                pos,
                                                pos_in_region,
                                                region_rc[pos], # read base
                                                region_seq[pos_in_region], # template base
                                                variant_read_base,
                                                variant_template_base,
                                                read_trinuc,
                                                template_trinuc] +
                                                baseinfo)
                            base_count_per_base_file.write(outline)
                            base_count_per_base_file.write("\n")



        ########################################################################
        ## EXACTLY X reads per base
        print(region_seq)
        print(len(region_seq))

        print("BC1exact")
        print(base_counter1_exact.shape)
        for pos, data in enumerate(base_counter1_exact):
            aml_loc = int(pos in r1_targets) + int(pos in r1_aml_only)
            if aml_loc==2:
                variant_read_base = var_base
                variant_template_base = var_base
            else:
                variant_read_base = '-'
                variant_template_base = '-'
            pos_in_region = pos

            if pos==0:
                read_trinuc='X'+region_seq[0:(pos+2)]
                template_trinuc='X'+region_seq[0:(pos_in_region+2)]
            elif pos==(len(region_seq)-1):
                read_trinuc=region_seq[(pos-1):(pos+2)]+'X'
                template_trinuc=region_seq[(pos_in_region-1):(pos_in_region+2)]+'X'
            else:
                read_trinuc=region_seq[(pos-1):(pos+2)]
                template_trinuc=region_seq[(pos_in_region-1):(pos_in_region+2)]

            print("Position: %d: read_trinuc: %s" % (pos,read_trinuc))
            print("Position: %d: template_trinuc: %s" % (pos,template_trinuc))

            for cindex,cov in enumerate(EXACT_COVERS):
                for nindex,nuc in enumerate(NUCS):
                    count = data[(cindex,nindex)]
                    total = np.sum(data[cindex])

                    alt_read_trinuc = read_trinuc[0] + nuc + read_trinuc[2]
                    alt_template_trinuc = template_trinuc[0] + nuc + template_trinuc[2]

                    baseinfo = [alt_read_trinuc,alt_template_trinuc,nuc,'.',total,count,'.',cov,'.'] 

                    outline = tabprint([aml_loc,
                                        region_index,
                                        refiter,
                                        region_strand,
                                        1,
                                        pos,
                                        pos_in_region,
                                        region_seq[pos], # read base
                                        region_seq[pos_in_region], # template base
                                        variant_read_base,
                                        variant_template_base,
                                        read_trinuc,
                                        template_trinuc] +
                                        baseinfo)
                    base_count_per_base_file.write(outline)
                    base_count_per_base_file.write("\n")

        print("BC2exact")
        print(base_counter2_exact.shape)
        for pos, data in enumerate(base_counter2_exact):
            aml_loc = int(pos in r2_targets) + int(pos in r2_aml_only)
            if aml_loc==2:
                variant_read_base = reverse_complement(var_base)
                variant_template_base = var_base
            else:
                variant_read_base = '-'
                variant_template_base = '-'
            pos_in_region = region_len - pos - 1

            if pos==0:
                read_trinuc='X'+region_rc[0:(pos+2)]
                template_trinuc=region_seq[(pos_in_region-1):(pos_in_region+2)]+'X'
            elif pos==(len(region_seq)-1):
                read_trinuc=region_rc[(pos-1):(pos+2)]+'X'
                template_trinuc='X'+region_seq[0:(pos_in_region+2)]
            else:
                read_trinuc=region_rc[(pos-1):(pos+2)]
                template_trinuc=region_seq[(pos_in_region-1):(pos_in_region+2)]
                print("Position: %d: read_trinuc: %s" % (pos,read_trinuc))
                print("Position: %d: template_trinuc: %s" % (pos,template_trinuc))

            for cindex,cov in enumerate(EXACT_COVERS):
                for nindex,nuc in enumerate(NUCS):
                    count = data[(cindex,nindex)]
                    total = np.sum(data[cindex])

                    alt_read_trinuc = read_trinuc[0] + nuc + read_trinuc[2]
                    alt_template_trinuc = template_trinuc[0] + reverse_complement(nuc) + template_trinuc[2]

                    baseinfo = [alt_read_trinuc,alt_template_trinuc,nuc,'.',total,count,'.',cov,'.'] 

                    outline = tabprint([aml_loc,
                                        region_index,
                                        refiter,
                                        region_strand,
                                        2,
                                        pos,
                                        pos_in_region,
                                        region_rc[pos], # read base
                                        region_seq[pos_in_region], # template base
                                        variant_read_base,
                                        variant_template_base,
                                        read_trinuc,
                                        template_trinuc] +
                                        baseinfo)
                    base_count_per_base_file.write(outline)
                    base_count_per_base_file.write("\n")

        ## AT LEAST X reads per base
        print("BC1over")
        print(base_counter1.shape)
        for pos, data in enumerate(base_counter1):
            aml_loc = int(pos in r1_targets) + int(pos in r1_aml_only)
            if aml_loc==2:
                variant_read_base = var_base
                variant_template_base = var_base
            else:
                variant_read_base = '-'
                variant_template_base = '-'
            pos_in_region = pos

            if pos==0:
                read_trinuc='X'+region_seq[0:(pos+2)]
                template_trinuc='X'+region_seq[0:(pos_in_region+2)]
            elif pos==(len(region_seq)-1):
                read_trinuc=region_seq[(pos-1):(pos+2)]+'X'
                template_trinuc=region_seq[(pos_in_region-1):(pos_in_region+2)]+'X'
            else:
                read_trinuc=region_seq[(pos-1):(pos+2)]
                template_trinuc=region_seq[(pos_in_region-1):(pos_in_region+2)]
            print("Position: %d: read_trinuc: %s" % (pos,read_trinuc))
            print("Position: %d: template_trinuc: %s" % (pos,template_trinuc))

            for cindex,cov in enumerate(MIN_COVERS):
                for nindex,nuc in enumerate(NUCS):
                    count = data[(cindex,nindex)]
                    total = np.sum(data[cindex])

                    alt_read_trinuc = read_trinuc[0] + nuc + read_trinuc[2]
                    alt_template_trinuc = template_trinuc[0] + nuc + template_trinuc[2]

                    baseinfo = [alt_read_trinuc,alt_template_trinuc,nuc,'.',total,count,cov,'.','.'] 

                    outline = tabprint([aml_loc,
                                        region_index,
                                        refiter,
                                        region_strand,
                                        1,
                                        pos,
                                        pos_in_region,
                                        region_seq[pos], # read base
                                        region_seq[pos_in_region], # template base
                                        variant_read_base,
                                        variant_template_base,
                                        read_trinuc,
                                        template_trinuc] +
                                        baseinfo)
                    base_count_per_base_file.write(outline)
                    base_count_per_base_file.write("\n")

        print("BC2over")
        print(base_counter2.shape)
        for pos, data in enumerate(base_counter2):
            aml_loc = int(pos in r2_targets) + int(pos in r2_aml_only)
            if aml_loc==2:
                variant_read_base = reverse_complement(var_base)
                variant_template_base = var_base
            else:
                variant_read_base = '-'
                variant_template_base = '-'
            pos_in_region = region_len - pos - 1

            if pos==0:
                read_trinuc='X'+region_rc[0:(pos+2)]
                template_trinuc=region_seq[(pos_in_region-1):(pos_in_region+2)]+'X'
            elif pos==(len(region_seq)-1):
                read_trinuc=region_rc[(pos-1):(pos+2)]+'X'
                template_trinuc='X'+region_seq[0:(pos_in_region+2)]
            else:
                read_trinuc=region_rc[(pos-1):(pos+2)]
                template_trinuc=region_seq[(pos_in_region-1):(pos_in_region+2)]
            print("Position: %d: read_trinuc: %s" % (pos,read_trinuc))
            print("Position: %d: template_trinuc: %s" % (pos,template_trinuc))

            for cindex,cov in enumerate(MIN_COVERS):
                for nindex,nuc in enumerate(NUCS):
                    count = data[(cindex,nindex)]
                    total = np.sum(data[cindex])

                    alt_read_trinuc = read_trinuc[0] + nuc + read_trinuc[2]
                    alt_template_trinuc = template_trinuc[0] + reverse_complement(nuc) + template_trinuc[2]

                    baseinfo = [alt_read_trinuc,alt_template_trinuc,nuc,'.',total,count,cov,'.','.'] 

                    outline = tabprint([aml_loc,
                                        region_index,
                                        refiter,
                                        region_strand,
                                        2,
                                        pos,
                                        pos_in_region,
                                        region_rc[pos], # read base
                                        region_seq[pos_in_region], # template base
                                        variant_read_base,
                                        variant_template_base,
                                        read_trinuc,
                                        template_trinuc] +
                                        baseinfo)
                    base_count_per_base_file.write(outline)
                    base_count_per_base_file.write("\n")



    ########################################################################
    # Write the within tag errors to pickle file and to table
    # Just save this for the first reference sequence.... for now
    logger.info('Writing within tag error to a file')
    WITHIN_TAG_ERRS = list_WITHIN_TAG_ERRS[0]
    pickle.dump(
        WITHIN_TAG_ERRS,
        open(args.within_tag_errors, 'wb'),
        pickle.HIGHEST_PROTOCOL,
    )
    for i,pt in enumerate(PAIRED_TRINUCS):
        withintagerrors_file.write(tabprint(['R1']+list(pt)+list(WITHIN_TAG_ERRS[i,:,0]))+'\n')
        withintagerrors_file.write(tabprint(['R2']+list(pt)+list(WITHIN_TAG_ERRS[i,:,1]))+'\n')

    #######################################################################
    # Close output files
    base_count_per_pos_file.close()
    base_count_per_base_file.close()
    alignment_rate_file.close()
    withintagerrors_file.close()
    unaligned_reads_file.close()
