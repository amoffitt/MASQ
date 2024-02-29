import numpy as np
import sys
import os
import gzip
import pickle
import subprocess
import time
import re
from collections import Counter, defaultdict
import yaml
import datetime
import pprint
import pysam

from loguru import logger

from masq.primer_design.enzymes import load_enzyme_descriptors, \
    process_enzyme_cut_sites

from primer_design_functions import *


###############################################################################
# Time entire script
start0=time.time()

###############################################################################
# Load config file as first command line argument
configfile=sys.argv[1]
config=yaml.load(open(configfile),Loader=yaml.SafeLoader)
pprint.pprint(config)
sys.stdout.flush()
###############################################################################
# Load reference genome
print('loading reference genome pickle')
seq_pickle = config['ref_pickle']
seq_dic = pickle.load(open(seq_pickle, 'rb'))

###############################################################################
# Enzyme files
enz_folder= config['folder_with_cut_site_files']
enzymes= config['enzyme_list']
genomebuild=config["genomebuildforcutsites"]
enzyme_pos_fns = [os.path.join(enz_folder,x+"."+genomebuild+".sort.gz") for x in enzymes]

###############################################################################
# Enzyme cut site offsets and recognition sites
# 1st column is enzyme name, 2nd column is motif, 3rd column is motif with cut, 4th column is cut offset
cut_site_file = config['cutsite_offset_file']
# cut_offsets=dict()
# motifs=dict()
# recogsites=dict()
# with open(cut_site_file,'r') as f:
#    for line in f:
#         e=line.strip().split()[0]
#         motifs[e]=line.strip().split()[1]
#         recogsites[e]=line.strip().split()[2]
#         cut_offsets[e]=int(line.strip().split()[3])

enzyme_descriptors = load_enzyme_descriptors(cut_site_file, enzymes)

###############################################################################
# Process enzyme cut sites into dictionary
# Top level keys: enzymes; Next level keys: chromosome
print("Loading enzyme cut site information")
start = time.time()
cut_site_pickle_top = config['enzyme_top_pickle']
cut_site_pickle_btm = config['enzyme_btm_pickle']
reload_enzymes = True
sys.stdout.flush()

# logger.info("try loading enzyme cut sites from pickle files")

# if os.path.exists(cut_site_pickle_top) and os.path.isfile(cut_site_pickle_top) \
#         and os.path.exists(cut_site_pickle_btm) \
#         and os.path.isfile(cut_site_pickle_btm):
#     cut_sites_top = pickle.load(open(cut_site_pickle_top, 'rb'))
#     cut_sites_btm = pickle.load(open(cut_site_pickle_btm, 'rb'))

#     for x in enzymes:  # check that all of the enzymes are present in the pickle file
#         if ( (x not in cut_sites_top) or (x not in cut_sites_btm) ):
#             reload_enzymes = True
#             break
# else:
#     reload_enzymes = True

# logger.info("pickle files loaded")

if reload_enzymes:
    cut_sites_top, cut_sites_btm = process_enzyme_cut_sites(
        [enzyme_descriptors[ename] for ename in enzymes],
        enz_folder, genomebuild)

    # cut_sites_top = {}
    # cut_sites_btm = {}
    # for x in enzymes:
    #     cut_sites_top[x]=dict()
    #     cut_sites_btm[x]=dict()
    # # Loop over enzyme files and store cut site information
    # c = 0
    # for ename, efn in zip(enzymes, enzyme_pos_fns):
    #     print(ename)
    #     print(efn)
    #     sys.stdout.flush()

    #     edesc = enzyme_descriptors[ename]
    #     off_top = edesc.cut_offset
    #     off_btm = len(edesc.motif) - off_top

    #     with gzip.open(efn, 'rt') as f:
    #         for line in f:
    #             entries = line.strip().split()
    #             chrom = entries[0]
    #             pos_top = int(entries[1]) + off_top
    #             pos_btm = int(entries[1]) + off_btm

    #             if chrom in cut_sites_top[ename].keys():
    #                 cut_sites_top[ename][chrom].append(pos_top)
    #                 cut_sites_btm[ename][chrom].append(pos_btm)
    #             else:
    #                 cut_sites_top[ename][chrom]=[pos_top]
    #                 cut_sites_btm[ename][chrom]=[pos_btm]
    #             c = c+1
    #             #if (c % 10000)==0:
    #             #    print(tabprint([chrom,pos_top]))

    # Save to python file
    # logger.info("saving enzyme cut sites to pickle files")
    # pickle.dump(
    #     cut_sites_top, open(cut_site_pickle_top, 'wb'),
    #     pickle.HIGHEST_PROTOCOL)
    # pickle.dump(
    #     cut_sites_btm, open(cut_site_pickle_btm, 'wb'),
    #     pickle.HIGHEST_PROTOCOL)
    # logger.info("done saving enzyme cut sites to pickle files")

end = time.time()
print("Time elapsed: %0.2f" % (end-start))
sys.stdout.flush()

###############################################################################
# Load SNPs into dictionary with initial information
# Chrom, pos, ref, alt, strand (optional: top or bottom)
# If strand is included - ref and alt are expected to be flipped for bottom
print("Loading SNP info"); start = time.time()
snp_file = config['variant_file']
snpdict=dict()

with open(snp_file,'rt') as f:
    for line in f:
        entries=line.strip().split()
        chrom=entries[0]
        pos=int(entries[1])
        ref=entries[2]
        alt=entries[3]
        try:
            strand=entries[4]
        except:
            strand=''

        # check reference base
        if strand=='bottom':
            seqref=reverseComplement(seq_dic[chrom][pos-1])
        else: # top or no strand provided
            seqref=seq_dic[chrom][pos-1]
        if seqref is not ref:
            print("Error: reference base (and strand) provided do not match reference seq dictionary")
            print("%s: %d, %s" % (chrom,pos,strand))
            sys.exit()

        if strand=='top':
            ref_trinuc=seq_dic[chrom][pos-2:pos+1]
            alt_trinuc=ref_trinuc[0]+alt+ref_trinuc[2]

            snpid='_'.join([chrom,str(pos),strand])
            snpdict[snpid]=dict()
            initial_snp_dict(snpdict,snpid,chrom,pos,strand,ref,alt,ref_trinuc,alt_trinuc)

        elif strand=='bottom':
            ref_trinuc=reverseComplement(seq_dic[chrom][pos-2:pos+1])
            alt_trinuc=ref_trinuc[0]+alt+ref_trinuc[2]

            snpid='_'.join([chrom,str(pos),strand])
            snpdict[snpid]=dict()
            initial_snp_dict(snpdict,snpid,chrom,pos,strand,ref,alt,ref_trinuc,alt_trinuc)
        else:
            ref_trinuc_fwd=seq_dic[chrom][pos-2:pos+1]
            alt_trinuc_fwd=ref_trinuc_fwd[0]+alt+ref_trinuc_fwd[2]

            ref_trinuc_rev=reverseComplement(seq_dic[chrom][pos-2:pos+1])
            alt_trinuc_rev=ref_trinuc_rev[0]+reverseComplement(alt)+ref_trinuc_rev[2]

            ref_rev=reverseComplement(alt)
            alt_rev=reverseComplement(alt)

            snpid='_'.join([chrom,str(pos),'top'])
            snpdict[snpid]=dict()
            strand='top'
            initial_snp_dict(snpdict,snpid,chrom,pos,strand,ref,alt,ref_trinuc_fwd,alt_trinuc_fwd)

            snpid='_'.join([chrom,str(pos),'bottom'])
            snpdict[snpid]=dict()
            strand='bottom'
            initial_snp_dict(snpdict,snpid,chrom,pos,strand,ref_rev,alt_rev,ref_trinuc_rev,alt_trinuc_rev)

print_snp_dict(snpdict,False)
end = time.time(); print("Time elapsed: %0.2f" % (end-start))
sys.stdout.flush()

###############################################################################
# Load trinucleotides with high error rates
if config['filter_trinucleotides']:
    high_error_trinucs = list()
    with open(config['trinucleotide_file']) as f:
        for line in f:
            reftrinuc = line.strip().split()[0]
            alttrinuc = line.strip().split()[1]
            high_error_trinucs.append((reftrinuc,alttrinuc))
###############################################################################
# Filter snps against trinucleotides
if config['filter_trinucleotides']:
    print("Filtering on Ref/Alt Trinucleotides")
    for s in snpdict:
        if snpdict[s]['status'] is 'pass':
            trinucs = (snpdict[s]['ref_trinuc'], snpdict[s]['alt_trinuc'])
            print(s)
            print(trinucs)
            if trinucs in high_error_trinucs:
                print("found in high error list")
                snpdict[s]['status']='drop'
                snpdict[s]['drop_reason']='high_error_trinuc'

###############################################################################
# Check mutation changes against enzyme cut sites to make list of bad enzymes per mutation
# Recognition sites only need checking in top strand - if they are there, they are in bottom too
print("Checking mutations for cut sites")
bad_enzyme_choices=defaultdict(list)
for s in snpdict:
    if snpdict[s]['status'] is 'pass':
        print(s)
        chrom = snpdict[s]['chrom']
        pos = snpdict[s]['pos']
        alt = snpdict[s]['alt']
        strand  = snpdict[s]['strand']
        ref_context=seq_dic[chrom][pos-7:pos+6]
        if strand=='top':
            alt_context=seq_dic[chrom][pos-7:pos-1]+alt+seq_dic[chrom][pos:pos+6]
        else:
            alt_context=seq_dic[chrom][pos-7:pos-1]+reverseComplement(alt)+seq_dic[chrom][pos:pos+6]
        print(ref_context)
        print(alt_context)
        for ename, edesc in enzyme_descriptors.items():
            ref_hit = check_sequence_for_cut_site(ref_context, edesc.motif)
            alt_hit = check_sequence_for_cut_site(alt_context, edesc.motif)
            if (alt_hit and not ref_hit):
                print("%s: mutation introduces cut site for %s, %s" % (s,ename, edesc))
                bad_enzyme_choices[s].append(ename)
            if (ref_hit and not alt_hit):
                print("%s: mutation removes cut site for %s, %s" % (s,ename, edesc))
                bad_enzyme_choices[s].append(ename)
print("Bad enzyme choices due to cut sites introduced by mutations")
print(bad_enzyme_choices)

###############################################################################
# Check snps in target region for introducing cut sites
# Optionally dropping snps with other snps in nearby cut sites
print("Checking target region for SNPs and cut sites")
bam=config['wgs_bam'] # Can be single string or list of strings ! (must have same reference)
bam_ref_dic = pickle.load(open(config['wgs_ref'], 'rb'))
for s in snpdict:
    if snpdict[s]['status'] is 'pass':
        target_pos = snpdict[s]['pos']
        chrom=snpdict[s]['chrom']
        strand  = snpdict[s]['strand']
        if strand=='top':
            fulltargetstring="%s:%d-%d" % (snpdict[s]['chrom'],max(target_pos+config['frag_end_range'][0],1),target_pos+config['good_cut_range'][1])
        else:
            fulltargetstring="%s:%d-%d" % (snpdict[s]['chrom'],max(target_pos-config['good_cut_range'][1],1),target_pos-config['frag_end_range'][0])
        print(fulltargetstring)
        [snp_positions,seq_positions,snp_alt_bases,region_ref_seq] = snps_or_indels_in_region(bam,fulltargetstring,bam_ref_dic,basequal_cutoff=config['basequal_cutoff'],vaf_cutoff=config['vaf_cutoff'],indelaf_cutoff=config['indelaf_cutoff'],var_count_cutoff=config['var_count_cutoff'],indel_count_cutoff=config['indel_count_cutoff'])
        non_target_snps = [x for x in snp_positions if x!=target_pos]

        print(s)
        print(bam)
        print(fulltargetstring)
        print(snp_positions)
        print(snp_alt_bases)

        # Save indel positions to check and warn later
        snpdict[s]['indel_positions'] = [x for x,y in zip(snp_positions,snp_alt_bases) if y=='I' ]
        if len(snpdict[s]['indel_positions'])>0:
            print("Indel positions: %s" % s)
            print(snpdict[s]['indel_positions'])

        if config['drop_snps_in_full_target']:
            if len(non_target_snps)>0:
                print("%s: found SNPs in nearby region" % s)
                snpdict[s]['status']='drop'
                snpdict[s]['drop_reason']='snps_in_target_region'
                continue

        if config['drop_indel_in_full_target']:
            if 'I' in snp_alt_bases:
                snpdict[s]['status']='drop'
                snpdict[s]['drop_reason']='indel_in_target_region'
                continue
        else:
            for x,y in zip(non_target_snps,snp_alt_bases):
                ref_context=seq_dic[chrom][x-7:x+6]
                alt_context=seq_dic[chrom][x-7:x-1]+y+seq_dic[chrom][x:x+6]
                if (alt_hit and not ref_hit):
                    print("%s: SNP in target region introduces cut site for %s, %s" % (s,e,m))
                    bad_enzyme_choices[s].append(e)
                if (ref_hit and not alt_hit):
                    print("%s: SNP in target region removes cut site for %s, %s" % (s,e,m))
                    bad_enzyme_choices[s].append(e)
print("Bad enzyme choices due to cut sites introduced by SNPs")
print(bad_enzyme_choices)

###############################################################################
# Current available SNPs from previous filtering
pass_snpids= [x for x in snpdict.keys() if snpdict[x]['status']=='pass']
snp_chrom=[snpdict[x]['chrom'] for x in pass_snpids]
snp_pos=[snpdict[x]['pos'] for x in pass_snpids]
snp_strand=[snpdict[x]['strand'] for x in pass_snpids]

print("SNP DICT BEFORE ENZYME SELECTION")
for s in snpdict:
    if snpdict[s]['status']=='drop':
        print(s)
        print(snpdict[s]['drop_reason'])

# Make dictionaries of good and bad cuts for each enzyme/snp pair
print("Idenfying good and bad cut sites for each enzyme-snp pair"); start = time.time()
sys.stdout.flush()
good_cuts=dict()
bad_cuts=dict()
fragend_cuts=dict()
for e in enzymes:
    good_cuts[e]=dict()
    bad_cuts[e]=dict()
    fragend_cuts[e]=dict()

# Keep track if SNP has any possible enzymes with good cuts / no bad cuts
possible_enzyme_match_found = dict()

for e in enzymes:
    print(e)
    sys.stdout.flush()
    for chrom,pos,strand,snpid in zip(snp_chrom,snp_pos,snp_strand,pass_snpids):
        print("Collecting cut sites. Enzyme: %s. SNP %s" % (e,snpid))
        if chrom in cut_sites_top[e].keys():
            x_top = np.array(cut_sites_top[e][chrom])
            x_btm = np.array(cut_sites_btm[e][chrom])

            if strand=='top':
                # Intervals are all closed intervals, not overlapping
                bad_range = ( pos + int(config['bad_cut_range'][0]) , pos + int(config['bad_cut_range'][1]) )
                good_range = ( pos + int(config['good_cut_range'][0]) , pos + int(config['good_cut_range'][1]) )
                fragend_range = ( pos + int(config['frag_end_range'][0]) , pos + int(config['frag_end_range'][1]) )

                print("Ranges:")
                print(good_range)
                print(bad_range)

                bad_list = x_top[((x_top>=bad_range[0]) & (x_top<=bad_range[1]))]
                good_list = x_top[((x_top>=good_range[0]) & (x_top<=good_range[1]))]
                fragend_list = x_btm[((x_btm>=fragend_range[0]) & (x_btm<=fragend_range[1]))] 

                print("Cuts in range:")
                print(good_list)
                print(bad_list)

            else:
                bad_range = ( pos - int(config['bad_cut_range'][1]) , pos - int(config['bad_cut_range'][0]) )
                good_range = ( pos - int(config['good_cut_range'][1]) , pos - int(config['good_cut_range'][0]) )
                fragend_range = ( pos - int(config['frag_end_range'][1]) , pos - int(config['frag_end_range'][0]) )

                print("Ranges:")
                print(good_range)
                print(bad_range)

                bad_list = x_btm[((x_btm>=bad_range[0]) & (x_btm<=bad_range[1]))] 
                good_list = x_btm[((x_btm>=good_range[0]) & (x_btm<=good_range[1]))]
                fragend_list = x_top[((x_top>=fragend_range[0]) & (x_top<=fragend_range[1]))]

                print("Cuts in range:")
                print(good_list)
                print(bad_list)

            good_cuts[e][snpid]=good_list
            bad_cuts[e][snpid]=bad_list
            fragend_cuts[e][snpid]=fragend_list

            if ( (len(good_list)>0) and (len(bad_list)==0) ):
                print("possible match found: %s" % snpid)
                possible_enzyme_match_found[snpid]=1


        else: # no cuts on that chromosome
            good_cuts[e][snpid]=np.array([])
            bad_cuts[e][snpid]=np.array([])
            fragend_cuts[e][snpid]=np.array([])


end = time.time(); print("Time elapsed: %0.2f" % (end-start))
sys.stdout.flush()

###############################################################################
# How many SNPs can we query with this enzyme list?
def what_snps_with_this_enzyme_list(enzs,good_cuts,bad_cuts,fragend_cuts,available_snps):
    goodsnps=dict()
    fragends=dict()
    removesnps=[]
    minfragsize = int(config['PRIMER_PRODUCT_SIZE_RANGE'][0])
    print(enzs)
    for e in enzs:
        for s,cuts in good_cuts[e].items():
            if s in available_snps:
                if len(cuts)>0:
                    # add snp to list if there exist good cuts
                    if s in goodsnps:
                        goodsnps[s]=np.append(goodsnps[s],cuts)
                    else:
                        goodsnps[s]=cuts

        for s,cuts in bad_cuts[e].items():
            if s in available_snps:
                if (len(cuts)>0):
                    # remove snp if bad cuts are found with this enzyme
                    removesnps.append(s)

        for s,cuts in fragend_cuts[e].items():
            if s in available_snps:
                if s in fragends:
                    fragends[s]=np.append(fragends[s],cuts)
                else:
                    fragends[s]=cuts

    # check fragment size
    for s,goodcuts in goodsnps.items():
        if ( (s in fragends) and (len(fragends[s])>0)): # else no fragment ending cuts
            # calculate fragment size
            if snpdict[s]['strand']=='top':
                fsize= min(goodcuts) - max(fragends[s])
            else:
                fsize= min(fragends[s]) - max(goodcuts)
            print("Fragment size: %s" % s)
            print(fsize)
            # if smaller than allowed size, drop snp
            if fsize < minfragsize:
                print("Dropping %s" % s)
                removesnps.append(s)

    # check bad_enzyme_choices
    for s in goodsnps:
        for x in bad_enzyme_choices[s]:
            if x in enzs:
                removesnps.append(s)

    # check distance between good cuts
    for s,cuts in goodsnps.items():
        if s in available_snps:
            if len(cuts)>1:
                if snpdict[s]['strand']=='top':
                    dist_to_next=(np.sort(cuts)-min(cuts))[1]
                else:
                    dist_to_next=(max(cuts)-np.sort(cuts))[-2]
                if dist_to_next < config['drop_if_cut_site_dist_lt_x_bases']:
                    removesnps.append(s)

    for s in list(set(removesnps)):
        if s in goodsnps:
            del goodsnps[s]

    return goodsnps

###############################################################################
# Select enzymes in greedy approach - the one that gives the most snps when added
# Stop when target snp number is reached or adding enzymes doesn't help
print("Selecting enzymes that maximize snp list"); start = time.time()
print(enzymes)
sys.stdout.flush()

snps_curr = dict() # all snps from all batches
batch_num=1
available_snps = [x for x in snpdict.keys() if snpdict[x]['status']=='pass']
too_small_batch = []
enz_remain = list(enzymes)
enzymes_for_batch = dict()

max_snp_count= config['max_snp_select']
target_batch_size = int(np.ceil( config['target_batch_size'] * 2 )) # to account for dropped snps later
min_batch_size = config['min_batch_size']

while ((len(snps_curr)<max_snp_count) & (len(enz_remain)>0)): # keep selecting new snps

    # New batch
    # Reset things as necessary
    enz_curr = []
    enz_remain = list(enzymes) # all enzymes to start each batch
    snps_curr_batch = dict()
    print("Batch number: %d" % batch_num)
    print("Number of enzymes to test: %d" % len(enz_remain))
    print("Number of availalbe SNPs: %d" % len(available_snps))

    while ((len(snps_curr_batch)<target_batch_size) & (len(enz_remain)>0)): # keep adding to batch
        print("Target batch size: %d" % target_batch_size)
        print("Current batch size: %d" % len(snps_curr_batch))

        nmax=0 # start over with number of snps when go over enz again
        print("Starting next iteration to find 1 enzyme to add")
        for i,e in enumerate(enz_remain):
            print("Iterations %d: Enzyme %s" % (i,e))
            # Enzyme list to test, add one enzyme at a time
            enz_test = list(enz_curr)
            enz_test.append(e)

            goodsnps = what_snps_with_this_enzyme_list(enz_test,good_cuts,bad_cuts,fragend_cuts,available_snps)
            print(goodsnps)
            # is this the best enzyme addition we've seen?
            n = len(goodsnps)
            print("By adding %s, we find %d SNPs" % (e,n))
            if n>nmax:
                nmax=n
                imax=i
                snpsmax=goodsnps
                print("Better than before")
        if nmax>0:
            print("Best enzyme index: %d" % imax)
            print(enz_remain)
            print(enz_remain[imax])
            print(nmax)

        # is this enzyme addition better than without it?
        if nmax>len(snps_curr_batch):

            chosen_enz = enz_remain[imax]
            print("Enzyme: %s" % chosen_enz)
            enz_curr.append(chosen_enz)
            print("Previous SNPs")
            if len(snps_curr_batch)>0:
                prev=set(snps_curr_batch.keys())
                print(prev)
            print("New SNPs")
            new=set(snpsmax.keys())
            print(new)
            if len(snps_curr_batch)>0:
                print("Lost from old set")
                print(prev.difference(new))
                print("Added in new set")
                print(new.difference(prev))
            snps_curr_batch = snpsmax
            enz_remain.remove(chosen_enz) # remove i
            print("#####################################")
        else:
            break
        sys.stdout.flush()

    # Done with batch - adjust snp lists accordingly
    print("Batch %d final selection" % batch_num)

    if len(snps_curr_batch)>min_batch_size:
        if len(snps_curr_batch)<target_batch_size: # between min and target
            snps_curr.update(snps_curr_batch) # add current batch to full list
            # remove snps from available list
            # also add batch info to snp dict
            for s in snps_curr_batch.keys():
                available_snps.remove(s)
                snpdict[s]['batch']=batch_num
            enzymes_for_batch[batch_num] = enz_curr
            batch_num += 1
            new_snps_curr_batch = snps_curr_batch
        else: # batch is too big, only keep max number
            new_snps_curr_batch = {k: snps_curr_batch[k] for k in list(snps_curr_batch)[:target_batch_size]}
            snps_curr.update(new_snps_curr_batch) # add current batch to full list
            # remove snps from available list
            # also add batch info to snp dict
            for s in new_snps_curr_batch.keys():
                available_snps.remove(s)
                snpdict[s]['batch']=batch_num
            enzymes_for_batch[batch_num] = enz_curr 
            batch_num += 1
    else:
        too_small_batch.append(snps_curr_batch.keys())
        break

    print(new_snps_curr_batch)
    print("All selected SNPs")
    print(snps_curr)
    print("Remaining SNPs")
    print(available_snps)
    print("################################")

print("All selected SNPs")
print(snps_curr)
print("Remaining SNPs: %d" % len(available_snps))
end = time.time(); print("Time elapsed: %0.2f" % (end-start))
sys.stdout.flush()

###############################################################################
# Given final snp list get enzymes assignments and cut distances
# Update SNP dictionary with pass, drop, reasons etc
print("Collecting information on final snp and enzyme list"); start = time.time()
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
                failedsnp_goodcuts='; '.join([e+':'+','.join([str(x) for x in good_cuts[e][s]]) for e in enz_curr])
                failedsnp_badcuts='; '.join([e+':'+','.join([str(x) for x in bad_cuts[e][s]]) for e in enz_curr])
                failedsnp_fragcuts='; '.join([e+':'+','.join([str(x) for x in fragend_cuts[e][s]]) for e in enz_curr])
                #TODO different enzymes enz-curr for different batches

                snpdict[s]['good_cuts']=failedsnp_goodcuts
                snpdict[s]['bad_cuts']=failedsnp_badcuts
                snpdict[s]['fragend_cuts']=failedsnp_fragcuts
                snpdict[s]['status']='drop'
                snpdict[s]['drop_reason']='enzyme_cut_compatibility'

print_snp_dict(snpdict,True)
end = time.time(); print("Time elapsed: %0.2f" % (end-start))
sys.stdout.flush()

###############################################################################
#Run PRIMER3
date=datetime.datetime.now().strftime('%Y-%m-%d.%H-%M')
primer3file=config['output_folder']+"primer3.input."+config['sample']+"."+date+".txt"
blatqueryfile=config['output_folder']+"blat_query.fa."+config['sample']+"."+date+".txt"
blatresultfile=config['output_folder']+"blat_results.out."+config['sample']+"."+date+".txt"
print(primer3file)
print(blatqueryfile)
print(blatresultfile)
# make output directory if it doesn't exist
os.makedirs(os.path.dirname(primer3file), exist_ok=True)

primer3results=dict()
#primer3 = config['primer3']
primer3 = "primer3_core" # if installed in environment or on path

print("Running primer3")
sys.stdout.flush()

with open(blatqueryfile,'w') as blatf:
    for snpid,snpinfo in snpdict.items():
        if snpinfo['status']=='pass':
            print("snp")
            write_primer3_input_file(
                primer3file,
                snpid,
                snpdict[snpid]['target_seq_for_primer_search'],
                snpdict[snpid]['strand'],
                snpdict[snpid]['dist_from_mut_to_upstream_cut'],
                config)
            print("")
            p=subprocess.Popen("%s %s" % (primer3,primer3file), shell=True, stdout=subprocess.PIPE)
            primer3out, err = p.communicate()
            print(primer3out.decode('ascii'))
            print("")
            sys.stdout.flush()

            # Store all the primer3 results
            primer3results[snpid]=dict()
            for line in primer3out.decode('ascii').split('\n'):
                if line.startswith('PRIMER'):
                    t,val=line.split("=")
                    primer3results[snpid][t]=val

            if "PRIMER_PAIR_NUM_RETURNED=0" in primer3out.decode('ascii'):
                snpdict[snpid]['status']='drop'
                snpdict[snpid]['drop_reason']='primer3_nonefound'

                snpdict[snpid]['left_primer_explanation']=primer3results[snpid]['PRIMER_LEFT_EXPLAIN']
                snpdict[snpid]['right_primer_explanation']=primer3results[snpid]['PRIMER_RIGHT_EXPLAIN']
            elif "PRIMER_ERROR" in primer3out.decode('ascii'):
                snpdict[snpid]['status']='drop'
                snpdict[snpid]['drop_reason']='primer3_error_seelog'

            else: # primer pairs found!
                for i in range(config['PRIMER_NUM_RETURN']):
                    t="PRIMER_LEFT_%d_SEQUENCE" % i
                    if t in primer3results[snpid].keys():
                        seq=primer3results[snpid][t]
                        blatf.write(">"+snpid+"_"+t+"\n")
                        blatf.write(seq+"\n")

                    t="PRIMER_RIGHT_%d_SEQUENCE" % i
                    if t in primer3results[snpid].keys():
                        seq=primer3results[snpid][t]
                        blatf.write(">"+snpid+"_"+t+"\n")
                        blatf.write(seq+"\n")

###############################################################################
# Run BLAT on all primer options
start=time.time()
print("Running BLAT on all primer options for all SNPs")
sys.stdout.flush()
run_blat(blatqueryfile,blatresultfile,config)
end = time.time(); print("Time elapsed: %0.2f" % (end-start))
sys.stdout.flush()

###############################################################################
# Process BLAT results
# Counter for number of hits per sequence
blat_hits=dict()

with open(blatresultfile,'r') as blatr:
    for line in blatr:
        s=line.split()[9]
        s_split=s.split("_")
        snpid="_".join(s_split[0:3])
        leftright=s_split[4]
        primerid=s_split[5]
        gaps=int(line.split()[6])
        # Only count entry if score is X away from length of primer(query)
        plen=int(line.split()[10])
        score=int(line.split()[0])
        if (score>=(plen - config['blat_num_mismatches'])):
            print("%s - %d - %d - %d" % (s,plen,score,gaps))
            if snpid in blat_hits:
                if leftright in blat_hits[snpid]:
                    blat_hits[snpid][leftright].update([primerid])
                else:
                    blat_hits[snpid][leftright]=Counter()
                    blat_hits[snpid][leftright].update([primerid])
            else:
                blat_hits[snpid]=dict()
                blat_hits[snpid][leftright]=Counter()
                blat_hits[snpid][leftright].update([primerid])


###############################################################################

# Look for valid pairs of primers based on blat hits
# Also check for SNPs in primer pairs

valid_primer_pairs=dict()
best_primer_pair=dict()

for snpid,counts in blat_hits.items():
    perfectfound=False
    perfectleft =[]
    perfectright =[]
    okleft =[]
    okright =[]

    print(snpid)
    print(counts)
    print("")

    if 'LEFT' in counts:
        for pid,ct in counts['LEFT'].items():
            if ct==config['blat_perfect_num_hits']:
                perfectleft.append(pid)
                print("Perfect left: %s - %s" % (snpid,pid))
            elif ct<config['blat_ok_num_hits']:
                okleft.append(pid)
                print("OK left: %s - %s" % (snpid,pid))
    if 'RIGHT' in counts:
        for pid,ct in counts['RIGHT'].items():
            if ct==config['blat_perfect_num_hits']:
                perfectright.append(pid)
                print("Perfect right: %s - %s" % (snpid,pid))
            elif ct<config['blat_ok_num_hits']:
                okright.append(pid)
                print("OK right: %s - %s" % (snpid,pid))

    # check for perfect pair
    perfectpairs=list(set(perfectleft).intersection(perfectright))
    if len(perfectpairs)>0:
        valid_primer_pairs[snpid]=perfectpairs
        ps=valid_primer_pairs[snpid]
        # check perfect pairs for snps in primers before skipping next step
        for p in ps:
            (primerstringL,primerstringR)=get_primer_coordinates(p,snpid,primer3results,snpdict)
            [snp_positionsL,seq_positions,snp_alt_bases,region_ref_seq] = snps_or_indels_in_region(bam,primerstringL,bam_ref_dic,basequal_cutoff=config['basequal_cutoff'],vaf_cutoff=config['vaf_cutoff'],indelaf_cutoff=config['indelaf_cutoff'],var_count_cutoff=config['var_count_cutoff'],indel_count_cutoff=config['indel_count_cutoff'])
            [snp_positionsR,seq_positions,snp_alt_bases,region_ref_seq] = snps_or_indels_in_region(bam,primerstringR,bam_ref_dic,basequal_cutoff=config['basequal_cutoff'],vaf_cutoff=config['vaf_cutoff'],indelaf_cutoff=config['indelaf_cutoff'],var_count_cutoff=config['var_count_cutoff'],indel_count_cutoff=config['indel_count_cutoff'])
            if (len(snp_positionsL)>0) or (len(snp_positionsR)>0):
                valid_primer_pairs[snpid].remove(p)
                print("Found SNP in primer pair: %s" % p)
                print("Left primer: %s" % primerstringL)
                print("Right primer: %s" % primerstringR)
        if len(valid_primer_pairs[snpid])>0:
            perfectfound=True

    if not perfectfound:
    # check for one perfect and one ok
        ok_perf_pairs=list(set(perfectleft).intersection(okright))
        ok_perf_pairs.extend(list(set(perfectright).intersection(okleft)))

        best_pairs=[]
        m=config['blat_perfect_num_hits']+config['blat_ok_num_hits'] # min hits combined across 2 primers
        for p in ok_perf_pairs:
            m_obs = blat_hits[snpid]['LEFT'][p] + blat_hits[snpid]['RIGHT'][p]
            if m_obs<m:
                m=m_obs
                best_pairs=[p]
            elif m_obs==m:
                best_pairs.append(p)

        valid_primer_pairs[snpid]=best_pairs

    print(valid_primer_pairs[snpid])
    print("Has valid primer pairs before SNP checking")
    # Further selection based on product size (larger is better)
    if len(valid_primer_pairs[snpid])>0:
        m=0
        ps=valid_primer_pairs[snpid]
        ps.sort(key=float) # sort so ties are broken by lowest number, which has best score from primer3
        print(ps)

        ps_no_snps = list(ps)
        # Check valid primer pairs for snps in primer, drop if SNP in primer
        for p in ps:
            (primerstringL,primerstringR)=get_primer_coordinates(p,snpid,primer3results,snpdict)
            [snp_positionsL,seq_positions,snp_alt_bases,region_ref_seq] = snps_or_indels_in_region(bam,primerstringL,bam_ref_dic,basequal_cutoff=config['basequal_cutoff'],vaf_cutoff=config['vaf_cutoff'],indelaf_cutoff=config['indelaf_cutoff'],var_count_cutoff=config['var_count_cutoff'],indel_count_cutoff=config['indel_count_cutoff'])
            [snp_positionsR,seq_positions,snp_alt_bases,region_ref_seq] = snps_or_indels_in_region(bam,primerstringR,bam_ref_dic,basequal_cutoff=config['basequal_cutoff'],vaf_cutoff=config['vaf_cutoff'],indelaf_cutoff=config['indelaf_cutoff'],var_count_cutoff=config['var_count_cutoff'],indel_count_cutoff=config['indel_count_cutoff'])
            if (len(snp_positionsL)>0) or (len(snp_positionsR)>0):
                print("Found SNP in primer pair: %s" % p)
                print("Left primer: %s" % primerstringL)
                print("Right primer: %s" % primerstringR)
                ps_no_snps.remove(p)

        # Still has valid primer options
        if len(ps_no_snps)>0:
            m=0
            for p in ps_no_snps:
                prodsize=int(primer3results[snpid]["PRIMER_PAIR_%s_PRODUCT_SIZE" % p])
                print("%d: %s" % (prodsize,p))
                if prodsize>m:
                    m=prodsize
                    bestprimer=p
                    print("bigger product")
            best_primer_pair[snpid]=bestprimer
            print(bestprimer)
    else:
        snpdict[snpid]['status']='drop'
        snpdict[snpid]['drop_reason']='blat_hits'

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
