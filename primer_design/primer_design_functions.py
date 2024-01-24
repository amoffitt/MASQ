#! /usr/local/bin/python
import sys, time
import pysam
import argparse
import csv
import numpy as np
from scipy import stats
import re
import pickle
import subprocess

from masq.utils.io import tabprint

def reverseComplement(seq):
    INDICT = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N',
              'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n'}
    return "".join([INDICT[base] for base in seq[::-1]])

def base2int(base):
    BASES = ["N", "A", "C", "G", "T"]
    BASE2INT = dict([x[::-1] for x in enumerate(BASES)])
    return BASE2INT[base]

def int2base(intbase):
    BASES = ["A", "C", "G", "T"]
    return BASES[intbase]

def split_region_string(regionstring):
    chrom = regionstring.split(':')[0]
    start = int(regionstring.split(':')[1].split('-')[0])
    end = int(regionstring.split(':')[1].split('-')[1])
    return [chrom,start,end]

def check_chr_in_seqdict(chrom,seq_dic):
    if chrom not in seq_dic:
        chrom = chrom[3:]
    return chrom

def check_sequence_for_cut_site(sequence,pattern):
    results = re.findall(pattern,sequence)
    if len(results)>0:
        return True # pattern is in seuqence
    else:
        return False

def snps_or_indels_in_region(bamfiles,regionstring,seq_dic,basequal_cutoff=28,vaf_cutoff=0.05,indelaf_cutoff=0.05,var_count_cutoff=2,indel_count_cutoff=2):
    # Coordinates should be 0-based: (VCF/IGV position)-1 = Python Position
    # Regions are defined as [start,stop), including start but not stop positions
    # bamfiles here can be one string for one bam, or a list of strings for multiple bams
    if isinstance(bamfiles,str):
        bamlist=[bamfiles] # convert to iterable list
    else:
        bamlist=bamfiles

    # # Extract coordinates
    print(regionstring)
    [chrom,start,end] = split_region_string(regionstring)
    print(chrom)
    chrom = check_chr_in_seqdict(chrom,seq_dic)
    print(chrom)
    # Get reference sequence
    region_ref_seq = seq_dic[chrom][(start-1):end]
    print(region_ref_seq)

    # Initialize
    snp_positions = np.array([])
    snp_alt_bases = []

    for bam in bamlist:
        # Load BAM file
        bamfile = pysam.Samfile(bam,"r")
        print(bamfile)
        # Pileup columns
        for pileupcolumn in bamfile.pileup(chrom, int(start)-1, int(end)+1, truncate = True, stepper = 'nofilter', max_depth=10000000): # stepper=all filters out pcr duplicates, set stepper=nofilter to not filter pcr duplicates
            # What position are we at?
            currpos = pileupcolumn.reference_pos
            refbaseint = base2int(seq_dic[chrom][currpos])

            basecounts = np.array([0,0,0,0,0]);
            totalcount = 0
            indelcount = 0
            has_something = False

            for pileupread in pileupcolumn.pileups:
                if pileupread.indel != 0:
                    indelcount+=1
                    totalcount+=1
                    continue;
                if ((pileupread.alignment.qual is None) or (pileupread.query_position is None)):
                    continue;
                basequal = ord(pileupread.alignment.qual[pileupread.query_position])-33
                if basequal<basequal_cutoff:
                    continue; # skip reads with quality less than filter value
                mapqual=pileupread.alignment.mapping_quality
                if mapqual<10:
                    continue; # skip reads with quality less than filter value
                # Extract base from read at this position
                totalcount+=1
                base=pileupread.alignment.seq[pileupread.query_position]
                baseint = base2int(base)
                basecounts[baseint] += 1

            # Drop N's
            basecounts = basecounts[1:]
            # Check if this pileup column has snp
            if np.sum(basecounts)==0:
                base_ratios = basecounts / 1.0
            else:
                base_ratios = basecounts / np.sum(basecounts)
            # Drop reference base from ratios
            ref_zero_base_ratios = base_ratios
            ref_zero_base_ratios[refbaseint-1] = 0
            # how many non-ref bases are here
            non_ref_base_count = np.sum(ref_zero_base_ratios> vaf_cutoff)

            if non_ref_base_count>0:
                ref_zero_base_counts = basecounts
                ref_zero_base_counts[refbaseint-1] = 0
                if max(ref_zero_base_counts)>var_count_cutoff:
                    altbase = int2base(np.argmax(ref_zero_base_ratios))
                    has_something = True
            if indelcount>indel_count_cutoff:
                if float(indelcount)/totalcount > indelaf_cutoff:
                    altbase = 'I' # indel
                    has_something = True
            if has_something:
                snp_positions = np.append(snp_positions,currpos+1) # put back into igv 1-based cooridinates
                snp_alt_bases.append(altbase)



    # Return 2 lists
    # Positions of alterations in region
    # Altered base/indel at each position
    seq_positions = [int(x) for x in (snp_positions - start)]
    snp_positions = [int(x) for x in snp_positions]

    return [snp_positions,seq_positions,snp_alt_bases,region_ref_seq]



def initial_snp_dict(snpdict,snpid,chrom,pos,strand,ref,alt,reftrinuc,alttrinuc):
    snpdict[snpid]['chrom']=chrom
    snpdict[snpid]['pos']=pos
    snpdict[snpid]['strand']=strand
    snpdict[snpid]['ref']=ref
    snpdict[snpid]['alt']=alt
    snpdict[snpid]['ref_trinuc']=reftrinuc
    snpdict[snpid]['alt_trinuc']=alttrinuc
    snpdict[snpid]['status']='pass'

def print_snp_dict(snpdict,passonly):
    for snpid,info in snpdict.items():
        if (snpdict[snpid]['status'])=='pass':
            print(snpid)
            for key,val in info.items():
                print("%s: %s" % (str(key),str(val)))
            print("\n######################\n")
        else:
            if not passonly:
                print(snpid)
                for key,val in info.items():
                    print("%s: %s" % (str(key),str(val)))
                print("\n######################\n")
    sys.stdout.flush()


def write_primer3_input_file(fn,snpid,templateseq,strand,dist,config):
    # Prepare PRIMER3 input file and run PRIMER3 on each SNP
    # SEQUENCE_TARGET should be position, followed by length!
    primer3text="""SEQUENCE_ID=%s
SEQUENCE_TEMPLATE=%s
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=%d
PRIMER_MIN_SIZE=%d
PRIMER_MAX_SIZE=%d
PRIMER_PRODUCT_SIZE_RANGE=%d-%d
PRIMER_PRODUCT_OPT_SIZE=%d
PRIMER_MIN_TM=%d
PRIMER_MAX_TM=%d
PRIMER_OPT_TM=%d
PRIMER_PAIR_MAX_DIFF_TM=%d
PRIMER_MIN_GC=%d
PRIMER_MAX_GC=%d
PRIMER_MAX_HAIRPIN_TH=%d
PRIMER_MAX_POLY_X=%d
PRIMER_NUM_RETURN=%d
PRIMER_TM_FORMULA=0
PRIMER_SALT_CORRECTIONS=0
PRIMER_EXPLAIN_FLAG=1
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s
=""" % (snpid, templateseq, config['PRIMER_OPT_SIZE'], config['PRIMER_MIN_SIZE'], config['PRIMER_MAX_SIZE'],
        config['PRIMER_PRODUCT_SIZE_RANGE'][0], config['PRIMER_PRODUCT_SIZE_RANGE'][1], config['PRIMER_PRODUCT_OPT_SIZE'],
        config['PRIMER_MIN_TM'], config['PRIMER_MAX_TM'], config['PRIMER_OPT_TM'], config['PRIMER_PAIR_MAX_DIFF_TM'],
        config['PRIMER_MIN_GC'], config['PRIMER_MAX_GC'], config['PRIMER_MAX_HAIRPIN_TH'], config['PRIMER_MAX_POLY_X'],
        config['PRIMER_NUM_RETURN'], config['primer3_thermo_param_folder'])

    if strand=='top':
        forcetext="SEQUENCE_FORCE_RIGHT_START=%d\n" % (len(templateseq)-1)
        targettext="SEQUENCE_TARGET=%d,%d\n" % ( len(templateseq)-dist-1,  2  )
    else:
        forcetext="SEQUENCE_FORCE_LEFT_START=0\n"
        targettext="SEQUENCE_TARGET=%d,%d\n" % (dist-1, 2)
    primer3text_plusforce = targettext + forcetext + primer3text

    with open(fn,'w') as f:
        f.write(primer3text_plusforce)
    print("writing file")
    print(snpid)

def run_blat(inputfile,outputfile,config,mode='primers'):
    ref=config['ref_fa']
    #blat = config['blat']
    blat = "blat" # if installed via conda or on path
    if mode=='primers':
        cmd="%s %s %s %s -tileSize=%d -stepSize=%d -minIdentity=%d -minScore=%d -maxIntron=%d -noHead" % (blat,ref,inputfile,outputfile,config['tileSize'],config['stepSize'],config['minIdentity'],config['minScore'],config['maxIntron'])
    else: # full length blat
        cmd="%s %s %s %s -tileSize=%d -stepSize=%d -minIdentity=%d -minScore=%d -maxIntron=%d -noHead" % (blat,ref,inputfile,outputfile,config['tileSize'],config['stepSize'],config['minIdentity_full'],config['minScore_full'],config['maxIntron_full'])
    p=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    out, err = p.communicate()


def get_primer_coordinates(primerid,snpid,primer3results,snpdict):
    chrom = snpdict[snpid]['chrom']
    targetpos = snpdict[snpid]['targetseq_pos1']
    left1 = int(primer3results[snpid]["PRIMER_LEFT_%s" % primerid].split(',')[0])+ targetpos + 1
    left2 = int(primer3results[snpid]["PRIMER_LEFT_%s" % primerid].split(',')[0])+ targetpos + 1 + int(primer3results[snpid]["PRIMER_LEFT_%s" % primerid].split(',')[1])
    right1 = int(primer3results[snpid]["PRIMER_RIGHT_%s" % primerid].split(',')[0])+ targetpos + 2 - int(primer3results[snpid]["PRIMER_RIGHT_%s" % primerid].split(',')[1])
    right2 = int(primer3results[snpid]["PRIMER_RIGHT_%s" % primerid].split(',')[0])+ targetpos + 2

    left_coord = "%s:%d-%d" % (chrom,left1,left2)
    right_coord = "%s:%d-%d" % (chrom,right1,right2)
    return (left_coord,right_coord)
