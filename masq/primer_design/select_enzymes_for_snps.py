import sys
import os
import time
import datetime
import argparse
import pprint
from typing import Optional

import yaml

from masq.primer_design.enzymes import load_enzyme_descriptors, \
    process_enzyme_cut_sites
from masq.primer_design.snp import process_snps, print_snp_dict, \
    filter_high_error_trinucleotides, check_snps_for_enzyme_cut_sites, \
    check_snps_in_target_region_for_cut_sites, \
    select_good_bad_cuts_for_enzyme_snp_pair, \
    greedy_select_enzimes, \
    filter_for_batch_size_duplication_and_no_primers, \
    update_snpdict_with_primer_info, store_snpdict_final, \
    update_snplist_with_enzyme_selection
from masq.primer_design.primer3_helpers import run_primer3
from masq.primer_design.blat_helpers import run_blat, process_blat_results, \
    find_valid_pairs_of_primers_based_on_blat_hits, run_full_blat_query
from masq.utils.reference_genome import ReferenceGenome


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Selects enzymes for SNPs and designs primers"
    )
    parser.add_argument(
        "config_file",
        help="primer design config file",
    )
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    if argv is None:
        argv = sys.argv[1:]

    args = parse_args(argv)

    ###########################################################################
    # Time entire script
    start0 = time.time()

    ###########################################################################
    # Load config file as first command line argument
    with open(args.config_file) as infile:
        config = yaml.load(infile, Loader=yaml.SafeLoader)

    pprint.pprint(config)
    sys.stdout.flush()
    ###########################################################################
    # Load reference genome
    print('loading reference genome')
    ref_genome = ReferenceGenome(config['ref_genome'])
    ref_genome.open()


    ###########################################################################
    # Enzyme files
    enz_folder = config['folder_with_cut_site_files']
    enzymes = config['enzyme_list']
    genomebuild = config["genomebuildforcutsites"]
    enzyme_pos_fns = [
        os.path.join(enz_folder, x + "." + genomebuild+".sort.gz")
        for x in enzymes]

    ###########################################################################
    # Enzyme cut site offsets and recognition sites
    # 1st column is enzyme name, 2nd column is motif, 3rd column is motif with
    # cut,
    # 4th column is cut offset
    cut_site_file = config['cutsite_offset_file']
    enzyme_descriptors = load_enzyme_descriptors(cut_site_file, enzymes)

    ###########################################################################
    # Process enzyme cut sites into dictionary
    # Top level keys: enzymes; Next level keys: chromosome
    print("Loading enzyme cut site information")
    start = time.time()
    sys.stdout.flush()

    cut_sites_top, cut_sites_btm = process_enzyme_cut_sites(
        [enzyme_descriptors[ename] for ename in enzymes],
        enz_folder, genomebuild)

    end = time.time()
    print(f"Time elapsed: {(end - start):0.2f}")
    sys.stdout.flush()

    ###########################################################################
    # Load SNPs into dictionary with initial information
    # Chrom, pos, ref, alt, strand (optional: top or bottom)
    # If strand is included - ref and alt are expected to be flipped for bottom
    print("Loading SNP info")
    start = time.time()
    snp_file = config['variant_file']
    snpdict = process_snps(snp_file, ref_genome)
    print_snp_dict(snpdict, False)

    end = time.time()
    print(f"Time elapsed: {(end - start):0.2f}")
    sys.stdout.flush()

    ###########################################################################
    # Load trinucleotides with high error rates
    if config['filter_trinucleotides']:
        snpdict = filter_high_error_trinucleotides(
            config['trinucleotide_file'], snpdict)

    ###########################################################################
    # Check mutation changes against enzyme cut sites to make list of bad
    # enzymes per mutation
    # Recognition sites only need checking in top strand - if they are there,
    # they are in bottom too
    print("Checking mutations for enzyme cut sites")
    bad_enzyme_choices = check_snps_for_enzyme_cut_sites(
        snpdict, enzyme_descriptors, ref_genome
    )

    ###########################################################################
    # Check snps in target region for introducing cut sites
    # Optionally dropping snps with other snps in nearby cut sites
    print("Checking target region for SNPs and cut sites")
    # Can be single string or list of strings ! (must have same reference)
    bam = config['wgs_bam']

    bad_enzyme_choices = check_snps_in_target_region_for_cut_sites(
        snpdict,
        enzyme_descriptors,
        ref_genome,
        bam,
        config,
        bad_enzyme_choices
    )
    print("Bad enzyme choices due to cut sites introduced by SNPs")
    print(bad_enzyme_choices)

    # #########################################################################
    # Current available SNPs from previous filtering
    # Make dictionaries of good and bad cuts for each enzyme/snp pair
    # in good_cuts, bad_cuts, fragend_cuts
    print("Idenfying good and bad cut sites for each enzyme-snp pair")
    start = time.time()
    sys.stdout.flush()

    # Keep track if SNP has any possible enzymes with good cuts / no bad cuts
    # in possible_enzyme_match_found

    good_cuts, bad_cuts, fragend_cuts, possible_enzyme_match_found = \
        select_good_bad_cuts_for_enzyme_snp_pair(
            snpdict, enzymes, cut_sites_top, cut_sites_btm, config)

    end = time.time()
    print(f"Time elapsed: {(end - start):0.2f}")
    sys.stdout.flush()


    ###########################################################################
    # Select enzymes in greedy approach - the one that gives the most snps when
    # added. Stop when target snp number is reached or adding enzymes doesn't
    # help
    print("Selecting enzymes that maximize snp list")
    start = time.time()
    print(enzymes)
    sys.stdout.flush()

    snps_curr, enzymes_for_batch, too_small_batch = greedy_select_enzimes(
        snpdict, enzymes,
        good_cuts, bad_cuts, fragend_cuts,
        bad_enzyme_choices, config
    )

    end = time.time()
    print(f"Time elapsed: {(end - start):0.2f}")
    sys.stdout.flush()

    ###########################################################################
    # Given final snp list get enzymes assignments and cut distances
    # Update SNP dictionary with pass, drop, reasons etc
    print("Collecting information on final snp and enzyme list")
    start = time.time()

    update_snplist_with_enzyme_selection(
        snpdict,
        snps_curr,
        enzymes_for_batch,
        good_cuts,
        bad_cuts,
        fragend_cuts,
        possible_enzyme_match_found,
        too_small_batch,
        enzyme_descriptors,
        ref_genome,
        config,
    )

    print_snp_dict(snpdict, True)

    end = time.time()
    print(f"Time elapsed: {(end - start):0.2f}")
    sys.stdout.flush()

    ###########################################################################
    # Run PRIMER3
    date = datetime.datetime.now().strftime('%Y-%m-%d.%H-%M')
    primer3file = config['output_folder'] + "primer3.input." + \
        config['sample'] + "." + date + ".txt"
    blatqueryfile = config['output_folder'] + "blat_query.fa." + \
        config['sample'] + "." + date + ".txt"
    blatresultfile = config['output_folder'] + "blat_results.out." + \
        config['sample'] + "." + date + ".txt"
    print(primer3file)
    print(blatqueryfile)
    print(blatresultfile)
    # make output directory if it doesn't exist
    os.makedirs(os.path.dirname(primer3file), exist_ok=True)

    # should be installed in environment or on path
    primer3 = config.get('primer3', "primer3_core")

    print("Running primer3")
    sys.stdout.flush()

    primer3results = run_primer3(
        snpdict,
        blatqueryfile,
        primer3file,
        primer3,
        config
    )

    ###########################################################################
    # Run BLAT on all primer options
    start = time.time()
    print("Running BLAT on all primer options for all SNPs")
    sys.stdout.flush()
    run_blat(blatqueryfile, blatresultfile, config)
    end = time.time()
    print(f"Time elapsed: {(end - start):0.2f}")
    sys.stdout.flush()

    ###########################################################################
    # Process BLAT results
    # Counter for number of hits per sequence

    blat_hits = process_blat_results(blatresultfile, config)

    ###########################################################################

    # Look for valid pairs of primers based on blat hits
    # Also check for SNPs in primer pairs

    best_primer_pair = find_valid_pairs_of_primers_based_on_blat_hits(
        snpdict,
        primer3results,
        blat_hits,
        ref_genome,
        config
    )

    ###########################################################################
    # Final filtering for batch size, duplicates, no primers found
    #
    snpdict = filter_for_batch_size_duplication_and_no_primers(
        snpdict,
        best_primer_pair,
        config
    )

    ###########################################################################
    # Loop over current pass snps
    # Update dict values for status based on primer results
    print("Getting final snp info")

    snpdict = update_snpdict_with_primer_info(
        snpdict,
        primer3results,
        blat_hits,
        best_primer_pair
    )

    sys.stdout.flush()

    ###########################################################################
    # Full length BLAT filter
    blatqueryfile = \
        config['output_folder'] + \
        "blat_query.full_length.fa." + \
        config['sample'] + \
        "." + date + ".txt"
    blatresultfile = \
        config['output_folder'] + \
        "blat_results.full_length.out." + \
        config['sample'] + \
        "." + date + ".txt"

    snpdict, _ = run_full_blat_query(
        blatqueryfile,
        blatresultfile,
        snpdict,
        ref_genome,
        config
    )
    sys.stdout.flush()

    ###########################################################################
    # Print final SNP results
    print_snp_dict(snpdict, False)
    print("")
    # Write to an output file
    # One row per SNP
    # date=datetime.datetime.now().strftime('%Y-%m-%d.%H-%M')
    outputfile = \
        config['output_folder'] + \
        "snp_primer_design_results." + \
        config['sample'] + "." + date + ".txt"

    store_snpdict_final(outputfile, snpdict)

    ###########################################################################
    # Final time
    end = time.time()
    print(f"Time elapsed for entire script: {(end - start0):0.2f}")
    sys.stdout.flush()
