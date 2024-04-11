import sys
import os
import pickle
import time
import datetime
import pprint

import yaml
import numpy as np

from masq.primer_design.enzymes import load_enzyme_descriptors, \
    process_enzyme_cut_sites
from masq.primer_design.snp import process_snps, print_snp_dict, \
    filter_high_error_trinucleotides, check_snps_for_enzyme_cut_sites, \
    check_snps_in_target_region_for_cut_sites, \
    select_good_bad_cuts_for_enzyme_snp_pair, \
    greedy_select_enzimes, \
    filter_for_batch_size_duplication_and_no_primers, \
    update_snpdict_with_primer_info, store_snpdict_final
from masq.primer_design.primer3_helpers import run_primer3
from masq.primer_design.blat_helpers import run_blat, process_blat_results, \
    find_valid_pairs_of_primers_based_on_blat_hits, run_full_blat_query
from masq.utils.reference_genome import ReferenceGenome


###############################################################################
# Time entire script
start0 = time.time()

###############################################################################
# Load config file as first command line argument
configfile = sys.argv[1]
config = yaml.load(open(configfile), Loader=yaml.SafeLoader)
pprint.pprint(config)
sys.stdout.flush()
###############################################################################
# Load reference genome
print('loading reference genome pickle')
seq_pickle = config['ref_pickle']
seq_dic = pickle.load(open(seq_pickle, 'rb'))
ref_genome = ReferenceGenome(config['ref_genome'])
ref_genome.open()


###############################################################################
# Enzyme files
enz_folder = config['folder_with_cut_site_files']
enzymes = config['enzyme_list']
genomebuild = config["genomebuildforcutsites"]
enzyme_pos_fns = [
    os.path.join(enz_folder, x + "." + genomebuild+".sort.gz")
    for x in enzymes]

###############################################################################
# Enzyme cut site offsets and recognition sites
# 1st column is enzyme name, 2nd column is motif, 3rd column is motif with cut,
# 4th column is cut offset
cut_site_file = config['cutsite_offset_file']
enzyme_descriptors = load_enzyme_descriptors(cut_site_file, enzymes)

###############################################################################
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

###############################################################################
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

###############################################################################
# Load trinucleotides with high error rates
if config['filter_trinucleotides']:
    snpdict = filter_high_error_trinucleotides(
        config['trinucleotide_file'], snpdict)

###############################################################################
# Check mutation changes against enzyme cut sites to make list of bad enzymes
# per mutation
# Recognition sites only need checking in top strand - if they are there,
# they are in bottom too
print("Checking mutations for enzyme cut sites")
bad_enzyme_choices = check_snps_for_enzyme_cut_sites(
    snpdict, enzyme_descriptors, ref_genome
)

###############################################################################
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

# ###############################################################################
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


###############################################################################
# Select enzymes in greedy approach - the one that gives the most snps when
# added. Stop when target snp number is reached or adding enzymes doesn't help
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

###############################################################################
# Given final snp list get enzymes assignments and cut distances
# Update SNP dictionary with pass, drop, reasons etc
print("Collecting information on final snp and enzyme list")
start = time.time()
for s in snpdict:
    if s in snps_curr:
        # Get enzyme list for this SNP's batch:
        enz_curr = enzymes_for_batch[snpdict[s]['batch']]

        # Find enzyme and cut site upstream
        sel_enz=''
        min_cut=snpdict[s]['pos']+1000
        max_cut=snpdict[s]['pos']-1000

        for e in enz_curr:
            cuts=good_cuts[e][s]
            if len(cuts)>0:
                if snpdict[s]['strand']=='top':
                    m=min(cuts)
                    if m<min_cut:
                        min_cut=m
                        sel_enz=e
                else:
                    m=max(cuts)
                    if m>max_cut:
                        max_cut=m
                        sel_enz=e

        snpdict[s]['enzyme']=sel_enz
        edesc = enzyme_descriptors[sel_enz]
        snpdict[s]['enzyme_recog_site']=edesc.recognition_site
        if snpdict[s]['strand']=='top':
            snpdict[s]['nearest_upstream_cut']=min_cut
        else:
            snpdict[s]['nearest_upstream_cut']=max_cut

        # Find closest cut site downstream (default 300)
        downstreamcut_min=snpdict[s]['pos'] + config['frag_end_range'][0]
        downstreamcut_max=snpdict[s]['pos'] - config['frag_end_range'][0]
        for e in enz_curr:
            cuts=fragend_cuts[e][s]
            if len(cuts)>0:
                if snpdict[s]['strand']=='top':
                    m=max(cuts)
                    if m>downstreamcut_min:
                        downstreamcut_min=m
                else:
                    m=min(cuts)
                    if m<downstreamcut_max:
                        downstreamcut_max=m

        if snpdict[s]['strand']=='top':
            snpdict[s]['nearest_downstream_cut']=downstreamcut_min
        else:
            snpdict[s]['nearest_downstream_cut']=downstreamcut_max

        # More information
        snpdict[s]['dist_from_mut_to_upstream_cut']=np.abs(snpdict[s]['nearest_upstream_cut']-snpdict[s]['pos'])
        pos1=min(snpdict[s]['nearest_upstream_cut'],snpdict[s]['nearest_downstream_cut'])
        pos2=max(snpdict[s]['nearest_upstream_cut'],snpdict[s]['nearest_downstream_cut'])
        if ( (pos2-pos1) < int(config['PRIMER_PRODUCT_SIZE_RANGE'][0])):
            snpdict[s]['status']='drop'
            snpdict[s]['drop_reason']='amplicon_too_short'
        else:
            snpdict[s]['target_seq_for_primer_search']=seq_dic[snpdict[s]['chrom']][pos1:pos2]
            snpdict[s]['targetseq_pos1']=pos1
            snpdict[s]['targetseq_pos2']=pos2
            snpdict[s]['targetseq_coordinates']="%s:%d-%d" % (snpdict[s]['chrom'],pos1+1,pos2)

            passedsnp_goodcuts='; '.join([e+':'+','.join([str(x) for x in good_cuts[e][s]]) for e in enz_curr])
            passedsnp_badcuts='; '.join([e+':'+','.join([str(x) for x in bad_cuts[e][s]]) for e in enz_curr])
            passedsnp_fragcuts='; '.join([e+':'+','.join([str(x) for x in fragend_cuts[e][s]]) for e in enz_curr])

            snpdict[s]['good_cuts']=passedsnp_goodcuts
            snpdict[s]['bad_cuts']=passedsnp_badcuts
            snpdict[s]['fragend_cuts']=passedsnp_fragcuts
    else:

        # Check if SNP had no chance - no enzyme that works with it alone
        if snpdict[s]['status']=='pass': # didn't fail earlier steps
            if s not in possible_enzyme_match_found.keys():
                snpdict[s]['status']='drop'
                snpdict[s]['drop_reason']='no_single_enzyme_matches'
            elif s in too_small_batch:
                snpdict[s]['status']='drop'
                snpdict[s]['drop_reason']='batch_too_small_from_enzyme_selection'
            else:
                # enz_curr = enzymes_for_batch[snpdict[s]['batch']]

                failedsnp_goodcuts='; '.join([e+':'+','.join([str(x) for x in good_cuts[e][s]]) for e in enz_curr])
                failedsnp_badcuts='; '.join([e+':'+','.join([str(x) for x in bad_cuts[e][s]]) for e in enz_curr])
                failedsnp_fragcuts='; '.join([e+':'+','.join([str(x) for x in fragend_cuts[e][s]]) for e in enz_curr])
                #TODO different enzymes enz-curr for different batches

                snpdict[s]['good_cuts']=failedsnp_goodcuts
                snpdict[s]['bad_cuts']=failedsnp_badcuts
                snpdict[s]['fragend_cuts']=failedsnp_fragcuts
                snpdict[s]['status']='drop'
                snpdict[s]['drop_reason']='enzyme_cut_compatibility'

print_snp_dict(snpdict, True)

# update_snplist_with_enzyme_selection(
#     snpdict,
#     snps_curr,
#     enzymes_for_batch,
#     good_cuts,
#     bad_cuts,
#     fragend_cuts,
#     possible_enzyme_match_found,
#     too_small_batch,
#     enzyme_descriptors,
#     ref_genome,
#     config,
# )

end = time.time()
print(f"Time elapsed: {(end - start):0.2f}")
sys.stdout.flush()

###############################################################################
#Run PRIMER3
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

#primer3 = config['primer3']
primer3 = "primer3_core"  # if installed in environment or on path

print("Running primer3")
sys.stdout.flush()

primer3results = run_primer3(
    snpdict,
    blatqueryfile,
    primer3file,
    primer3,
    config
)

###############################################################################
# Run BLAT on all primer options
start = time.time()
print("Running BLAT on all primer options for all SNPs")
sys.stdout.flush()
run_blat(blatqueryfile, blatresultfile, config)
end = time.time()
print(f"Time elapsed: {(end - start):0.2f}")
sys.stdout.flush()

###############################################################################
# Process BLAT results
# Counter for number of hits per sequence

blat_hits = process_blat_results(blatresultfile, config)

###############################################################################

# Look for valid pairs of primers based on blat hits
# Also check for SNPs in primer pairs

best_primer_pair = find_valid_pairs_of_primers_based_on_blat_hits(
    snpdict,
    primer3results,
    blat_hits,
    ref_genome,
    config
)

###############################################################################
# Final filtering for batch size, duplicates, no primers found
#
snpdict = filter_for_batch_size_duplication_and_no_primers(
    snpdict,
    best_primer_pair,
    config
)

###############################################################################
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

###############################################################################
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

blat_hits = run_full_blat_query(
    blatqueryfile,
    blatresultfile,
    snpdict,
    ref_genome,
    config
)
sys.stdout.flush()

###############################################################################
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

###############################################################################
# Final time
end = time.time()
print(f"Time elapsed for entire script: {(end - start0):0.2f}")
sys.stdout.flush()
