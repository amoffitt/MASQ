import numpy as np
import sys
import os
import pickle
import subprocess
import time
from collections import Counter
import yaml
import datetime
import pprint
import copy

from masq.primer_design.enzymes import load_enzyme_descriptors, \
    process_enzyme_cut_sites
from masq.primer_design.snp import process_snps, print_snp_dict, \
    filter_high_error_trinucleotides, check_snps_for_enzyme_cut_sites, \
    snps_or_indels_in_region, check_snps_in_target_region_for_cut_sites, \
    select_good_bad_cuts_for_enzyme_snp_pair, \
    greedy_select_enzimes, \
    update_snplist_with_enzyme_selection
from masq.primer_design.primer3_helpers import run_primer3
from masq.primer_design.blat_helpers import run_blat, process_blat_results, \
    find_valid_pairs_of_primers_based_on_blat_hits
from masq.utils.reference_genome import ReferenceGenome

from primer_design_functions import *


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
print("Final filtering for batch size, duplicates, no primers found")
# Get batch size
batch_size_counter = Counter()
for snpid,snpinfo in snpdict.items():
    print(snpid)
    if snpinfo['status']=='pass':
        print("pass")
        if snpid in best_primer_pair:
            print("has primer pair")

            # Update batch size
            b = snpdict[snpid]['batch']
            batch_size_counter.update([b])
            # If we've reached max for this batch, drop remainder
            if batch_size_counter[b] > config['max_batch_size']:
                snpdict[snpid]['status']='drop'
                snpdict[snpid]['drop_reason']='max_batch_size_reached'
        # Drop anything without valid primer pair
        else:
            snpdict[snpid]['status']='drop'
            snpdict[snpid]['drop_reason']='SNP_in_primer'

# Find batches that are too small
dropbatch_small=[]
for b,batchsize in batch_size_counter.items():
    if batchsize < config['min_batch_size']:
        dropbatch_small.append(b)

# Drop small batches
for snpid,snpinfo in snpdict.items():
    print("SNP")
    print(snpid)
    if snpinfo['status']=='pass':
        b = snpdict[snpid]['batch']
        if b in dropbatch_small:
            snpdict[snpid]['status']='drop'
            snpdict[snpid]['drop_reason']='batch_size_too_small'

# Check for duplicates
passed_snps = [x for x in snpdict.keys() if snpdict[x]['status']=='pass']
for snpid,snpinfo in snpdict.items():
    print(snpid)
    if snpinfo['status']=='pass':
        # Check for same position on both strands in final list
        if snpdict[snpid]['strand']!=config['strand_preference']:
            opposite_snpid='_'.join(snpid.split('_')[0:2]+['top'])
            if opposite_snpid in passed_snps: # top strand for same snp is in final list, drop this one
                snpdict[snpid]['status']='drop'
                snpdict[snpid]['drop_reason']='other_strand_same_snp_in_final_list'

#

###############################################################################
# Loop over current pass snps
# Update dict values for status based on primer results
print("Getting final snp info")

for snpid,snpinfo in snpdict.items():
    if snpinfo['status']=='pass':
        # Get assigned primer id here
        primerid = best_primer_pair[snpid]
        print("Primer selected - %s: %s" % (snpid,primerid))
        # Access dictionary of primer3 results to get relevant info
        snpdict[snpid]['primerID']=primerid
        if snpdict[snpid]['strand']=='top':
            print("top strand")
            # Right primer is cut-adjacent
            snpdict[snpid]['cutadj_primerseq']=primer3results[snpid]["PRIMER_RIGHT_%s_SEQUENCE" % primerid]
            snpdict[snpid]['downstream_primerseq']=primer3results[snpid]["PRIMER_LEFT_%s_SEQUENCE" % primerid]
            snpdict[snpid]['cutadj_melting_temp']=primer3results[snpid]["PRIMER_RIGHT_%s_TM" % primerid]
            snpdict[snpid]['downstream_melting_temp']=primer3results[snpid]["PRIMER_LEFT_%s_TM" % primerid]
            snpdict[snpid]['cutadj_primer_length']=len(primer3results[snpid]["PRIMER_RIGHT_%s_SEQUENCE" % primerid])
            snpdict[snpid]['downstream_primer_length']=len(primer3results[snpid]["PRIMER_LEFT_%s_SEQUENCE" % primerid])
            snpdict[snpid]['amplicon_length']=primer3results[snpid]["PRIMER_PAIR_%s_PRODUCT_SIZE" % primerid]
            snpdict[snpid]['cutadj_blat_unique']=blat_hits[snpid]['RIGHT'][primerid]
            snpdict[snpid]['downstream_blat_unique']=blat_hits[snpid]['LEFT'][primerid]
            snpdict[snpid]['cutadj_primer_GC']=primer3results[snpid]["PRIMER_RIGHT_%s_GC_PERCENT" % primerid]
            snpdict[snpid]['downstream_primer_GC']=primer3results[snpid]["PRIMER_LEFT_%s_GC_PERCENT" % primerid]
            snpdict[snpid]['cutadj_primer_coordinates']="%s:%d-%d" % (snpdict[snpid]['chrom'],
                int(primer3results[snpid]["PRIMER_RIGHT_%s" % primerid].split(',')[0])+ snpdict[snpid]['targetseq_pos1'] - int(primer3results[snpid]["PRIMER_RIGHT_%s" % primerid].split(',')[1]) + 2,
                int(primer3results[snpid]["PRIMER_RIGHT_%s" % primerid].split(',')[0]) + snpdict[snpid]['targetseq_pos1'] + 1)
            snpdict[snpid]['downstream_primer_coordinates']="%s:%d-%d" % (snpdict[snpid]['chrom'],
                int(primer3results[snpid]["PRIMER_LEFT_%s" % primerid].split(',')[0])+ snpdict[snpid]['targetseq_pos1'] + 1,
                int(primer3results[snpid]["PRIMER_LEFT_%s" % primerid].split(',')[0])+int(primer3results[snpid]["PRIMER_LEFT_%s" % primerid].split(',')[1])+ snpdict[snpid]['targetseq_pos1'] )
        else:
            print("bottom strand")
            # Left primer is cut-adjacent
            snpdict[snpid]['cutadj_primerseq']=primer3results[snpid]["PRIMER_LEFT_%s_SEQUENCE" % primerid]
            snpdict[snpid]['downstream_primerseq']=primer3results[snpid]["PRIMER_RIGHT_%s_SEQUENCE" % primerid]
            snpdict[snpid]['cutadj_melting_temp']=primer3results[snpid]["PRIMER_LEFT_%s_TM" % primerid]
            snpdict[snpid]['downstream_melting_temp']=primer3results[snpid]["PRIMER_RIGHT_%s_TM" % primerid]
            snpdict[snpid]['cutadj_primer_length']=len(primer3results[snpid]["PRIMER_LEFT_%s_SEQUENCE" % primerid])
            snpdict[snpid]['downstream_primer_length']=len(primer3results[snpid]["PRIMER_RIGHT_%s_SEQUENCE" % primerid])
            snpdict[snpid]['amplicon_length']=primer3results[snpid]["PRIMER_PAIR_%s_PRODUCT_SIZE" % primerid]
            snpdict[snpid]['cutadj_blat_unique']=blat_hits[snpid]['LEFT'][primerid]
            snpdict[snpid]['downstream_blat_unique']=blat_hits[snpid]['RIGHT'][primerid]
            snpdict[snpid]['cutadj_primer_GC']=primer3results[snpid]["PRIMER_LEFT_%s_GC_PERCENT" % primerid]
            snpdict[snpid]['downstream_primer_GC']=primer3results[snpid]["PRIMER_RIGHT_%s_GC_PERCENT" % primerid]
            snpdict[snpid]['cutadj_primer_coordinates']="%s:%d-%d" % (snpdict[snpid]['chrom'],
                int(primer3results[snpid]["PRIMER_LEFT_%s" % primerid].split(',')[0])+ snpdict[snpid]['targetseq_pos1'] + 1,
                int(primer3results[snpid]["PRIMER_LEFT_%s" % primerid].split(',')[0])+ int(primer3results[snpid]["PRIMER_LEFT_%s" % primerid].split(',')[1])+ snpdict[snpid]['targetseq_pos1'] )
            snpdict[snpid]['downstream_primer_coordinates']="%s:%d-%d" % (snpdict[snpid]['chrom'],
                int(primer3results[snpid]["PRIMER_RIGHT_%s" % primerid].split(',')[0])+ snpdict[snpid]['targetseq_pos1'] - int(primer3results[snpid]["PRIMER_RIGHT_%s" % primerid].split(',')[1]) + 2,
                int(primer3results[snpid]["PRIMER_RIGHT_%s" % primerid].split(',')[0])+ snpdict[snpid]['targetseq_pos1'] + 1)

        # Get full amplicon region coordinates
        posa=min(int(snpdict[snpid]['cutadj_primer_coordinates'].split(':')[1].split('-')[0]),
                 int(snpdict[snpid]['cutadj_primer_coordinates'].split(':')[1].split('-')[1]),
                 int(snpdict[snpid]['downstream_primer_coordinates'].split(':')[1].split('-')[0]),
                 int(snpdict[snpid]['downstream_primer_coordinates'].split(':')[1].split('-')[1]))
        posb=max(int(snpdict[snpid]['cutadj_primer_coordinates'].split(':')[1].split('-')[0]),
                 int(snpdict[snpid]['cutadj_primer_coordinates'].split(':')[1].split('-')[1]),
                 int(snpdict[snpid]['downstream_primer_coordinates'].split(':')[1].split('-')[0]),
                 int(snpdict[snpid]['downstream_primer_coordinates'].split(':')[1].split('-')[1]))

        # check for indels in amplicon and add as warning:
        indels = snpdict[snpid]['indel_positions']
        for i in indels:
            if (int(i)>=int(posa) and int(i)<=int(posb)):
                snpdict[snpid]['warnings'] = 'indel_in_amplicon_check_cut_sites'

        # Add region to output
        snpdict[snpid]['full_region_coordinates']="%s:%s-%s" % (snpdict[snpid]['chrom'],str(posa),str(posb))

sys.stdout.flush()

###############################################################################
# Full length BLAT filter
blatqueryfile=config['output_folder']+"blat_query.full_length.fa."+config['sample']+"."+date+".txt"
blatresultfile=config['output_folder']+"blat_results.full_length.out."+config['sample']+"."+date+".txt"
with open(blatqueryfile,'w') as blatf:
    passed_snps = [x for x in snpdict.keys() if snpdict[x]['status']=='pass']
    for snpid,snpinfo in snpdict.items():
        print(snpid)
        if snpinfo['status']=='pass':
            full_region=snpdict[snpid]['full_region_coordinates']
            # GET FULL LENGTH SEQ
            full_region_seq = seq_dic[full_region.split(':')[0]][int(full_region.split(':')[1].split('-')[0]):int(full_region.split(':')[1].split('-')[1])]
            # save seq to dictionary
            snpdict[snpid]['full_region_seq']=full_region_seq
            # Write to blat query file
            blatf.write(">"+snpid+"\n")
            blatf.write(full_region_seq+"\n")

# Run full length blat on all sequences at once
start=time.time()
print("Running BLAT on full length sequences for all SNPs")
sys.stdout.flush()
run_blat(blatqueryfile,blatresultfile,config,'full_length')
end = time.time(); print("Time elapsed: %0.2f" % (end-start))
sys.stdout.flush()

# Process blat results
blat_hits=Counter()
with open(blatresultfile,'r') as blatr:
    for line in blatr:
        snpid=line.split()[9]

        gaps=int(line.split()[6])
        plen=int(line.split()[10])
        score=int(line.split()[0])

        # only count those that pass min score
        if score>=config['minScore_full']:
            blat_hits.update([snpid])
            if blat_hits[snpid]>config['blat_full_num_hits']:
                snpdict[snpid]['status']='drop'
                snpdict[snpid]['drop_reason']='full_len_blat'
print("Done with full length blat")
sys.stdout.flush()

###############################################################################
# Print final SNP results
print_snp_dict(snpdict,False)
print("")
# Write to an output file
# One row per SNP
# date=datetime.datetime.now().strftime('%Y-%m-%d.%H-%M')
outputfile=config['output_folder']+"snp_primer_design_results."+config['sample']+"."+date+".txt"

# Manually set output order
header=['status','drop_reason','batch','full_region_coordinates','full_region_seq','chrom','pos','strand','ref_trinuc','alt_trinuc','enzyme','enzyme_recog_site',
    'amplicon_length','dist_from_mut_to_upstream_cut','nearest_downstream_cut','nearest_upstream_cut',
    'good_cuts','bad_cuts','fragend_cuts','targetseq_coordinates','target_seq_for_primer_search','left_primer_explanation',
    'right_primer_explanation','primerID','cutadj_primer_coordinates','cutadj_primerseq','cutadj_primer_length','cutadj_melting_temp',
    'cutadj_primer_GC','cutadj_blat_unique','downstream_primer_coordinates','downstream_primerseq','downstream_primer_length','downstream_melting_temp',
    'downstream_primer_GC','downstream_blat_unique','warnings']

with open(outputfile,'w') as f:
    f.write(tabprint(['SNP_ID']+header)+"\n")
    for snpid,snpinfo in snpdict.items():
        f.write(snpid+"\t")
        for h in header:
            if h in snpinfo:
                f.write(str(snpinfo[h])+"\t")
            else:
                f.write("."+"\t")
        f.write("\n")

###############################################################################
# Final time
end = time.time(); print("Time elapsed for entire script: %0.2f" % (end-start0))
sys.stdout.flush()
