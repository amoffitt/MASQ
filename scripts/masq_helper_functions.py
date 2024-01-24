'''
Helper functions for MASQ pipeline
'''

import sys
import operator
import re
import logging

from typing import Any, TextIO
from collections import Counter, defaultdict

import numpy as np

import matplotlib.pyplot as plt

import editdistance


# Logging
def setup_logger(logfile,name):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(str(logfile))
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    return(logger)
    # debug, info, warning, error


def reverseComplement(seq: str) -> str:
    """Get reverse complement of sequence."""
    indict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
              'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return "".join([indict[base] for base in seq[::-1]])


## get complement of sequence
def complement(seq):
    INDICT = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N',
              'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n'}
    return "".join([INDICT[base] for base in seq])

# splits cigar string into tuples of (letter, number)
def convert_cigar_string(cigar):
    # splits cigar string into tuples of (letter, number)
    match = re.findall(r'(\d+)(\w)', cigar)
    return [(y, int(x)) for x, y in match]


def hamming(
    str1: str, str2: str, use_edit_distance: bool = False
) -> int:
    """Compute hamming distance."""
    # only compares up until end of shortest string
    ne = operator.ne
    if use_edit_distance:
        min_len = min(len(str1), len(str2))
        return int(editdistance.eval(str1[:min_len],str2[:min_len]))

    return int(sum(map(ne, str1, str2)))

## simple class for counting pairs
def double_counter():
    return [0, 0]

## apply hamming distance to pair of reads
def test_pair(a1, a2, b1, b2, use_edit_distance=False):
    h1 = hamming(a1, b1, use_edit_distance=use_edit_distance)
    h2 = hamming(a2, b2, use_edit_distance=use_edit_distance)
    return h1 + h2

def hamming_withNs(str1,str2,maxNratio=0.2,use_edit_distance=False):
    ne = operator.ne
    c = operator.countOf
    L = min(len(str1),len(str2))
    ncount = max(c(str1[:L],'N'),c(str2[:L],'N'))
    if use_edit_distance:
        h = editdistance.eval(str1[:L],str2[:L])
    else:
        h = sum(map(ne, str1, str2))
    if float(ncount)/L > maxNratio:
        return L
    else:
        return (h - ncount)

def test_pair_withNs(a1, a2, b1, b2, maxNratio=0.2,use_edit_distance=False):
    h1 = hamming_withNs(a1, b1, use_edit_distance=use_edit_distance)
    h2 = hamming_withNs(a2, b2, use_edit_distance=use_edit_distance)
    return h1 + h2


########################################################################

# '''
#    Class for walking through a paired-end data of fastq.gz files
#    where the reads are organized in pairs in the file.
#    Once initialized with two filenames (for read1 and read2 fastq.gz files)
#    it returns an iterator which returns one read pair at a time.
# '''
# class readWalker():
#     ## standard  function for initilizing a class
#     def __init__(self, read1filename, read2filename):
#         self.read1file = gzip.open(read1filename, 'rt') # changed from r to rt
#         self.read2file = gzip.open(read2filename, 'rt')
#         gzip.open
#     ## standard function for classes that implement iteration
#     def __iter__(self):
#         return self

#     ## return the next pair
#     # Python2->3 replace next with __next__
#     # # Python2->3 replace next with readline
#     def __next__(self):
#         ## try to load data from read files/
#         ## fail at the end of file, so stop iterating at "except"
#         try:
#             ## load the read name (redundantly) from each file
#             pairID = self.read1file.readline().strip()
#             pairID = self.read2file.readline().strip()
#             ## read the sequence data
#             seq1            = self.read1file.readline().strip()
#             seq2            = self.read2file.readline().strip()
#             ## skip the "+" line
#             self.read1file.readline()
#             self.read2file.readline()
#             ## read the qualities
#             qual1           = self.read1file.readline().strip()
#             qual2           = self.read2file.readline().strip()
#             ## build the output object
#             nextPair        = [pairID,
#                                seq1, qual1,
#                                seq2, qual2]
#         ## if we failed to read, set nextPair to "None"
#             if not pairID:  # returned empty string
#                 nextPair = None
#         except:
#             nextPair = None

#         ## if we had to stop, close files and raise StopIteration
#         if nextPair == None:
#             self.close()
#             #return  # tried this to get it to stop
#             raise StopIteration
#         ## otherwise return the pair
#         else:
#             return nextPair
#     ## close the iterator
#     def close(self):
#         self.read1file.close()
#         self.read2file.close()

########################################################################


# takes as input taglist
# returns match dictionary
# needs further processing to update any previous counters
def cluster_rollup1(tag_list, show_progress=False):
    ## tag_list is an unsorted list of observed tags
    ## first statement gets unique elements, their count and a map from tag_list to unique_list
    unique_list, unique_inverse, unique_count = np.unique(tag_list, return_inverse=True, return_counts=True)
    ## assumes all tags are of equal length
    TAGLEN = len(unique_list[0])
    ## determine the order based on count, largest first
    ## NOTE: ADD JITTER HERE?
    order = np.argsort(-unique_count)
    order_reverse = np.argsort(order)
    ## reorder the unique elements to reflect count order
    unique_inverse = order_reverse[unique_inverse]
    unique_list    = unique_list[order]
    unique_count   = unique_count[order]
    L = len(unique_list)
    ## match_dict stores the first hamming-radius match for each entry (initialized to itself)
    match_dict = np.arange(L)
    ## omit dictionary is a dictionary for each position 0 <= x < TAGLEN
    ## such that when index == index,
    ## omit_dictionary[i] has keys of length TAGLEN -1
    ## derived from the previous entries tag[k] for k < index
    ## by omitting the i-th letter,
    ## and the value is the first index with that key.
    omit_dictionary = [{} for _ in range(TAGLEN)]
    counter = 0
    for index, tag in enumerate(unique_list):
        # print(tag)
        counter += 1
        if show_progress and counter % 100000 == 0:
            print(counter)
        ## keep a set of values where the i-th omit matches a previous key.
        match       = set()
        for i in range(TAGLEN):
            # print(i)
            newtag = tag[:i] + tag[(i+1):]
            try:
                ## check if it's already in the dictionary
                has_match = omit_dictionary[i][newtag]
                ## if it is, add it to the list of matches
                match.add(has_match)
            except KeyError:
                ## if not, there is a error raised
                ## and we now assign the present index to this entry
                omit_dictionary[i][newtag] = index
        ## THIS IS AT PRESENT A "STRICT" ROLLUP
        ## You are assigned to your top match by a sequence of hamming-1 distances
        ## that connect lower counts to higher counts
        matches = list(np.sort(list(match))) + [index]
        match_dict[index]  = match_dict[matches[0]] # select lowest index matching tag - then go look up what value it has stored in the dictionary
        # pprint.pprint(match_dict)
    return unique_list, unique_count, match_dict, match_dict[unique_inverse]


########################################################################
# takes as input a:
# 1) counter with tags and tag counts
# 2) dictionary with tags as keys, and counter as values ((read1,read2) as key, count as value)
# returns updated versions of these

def cluster_rollup2(vt_counter,vt_seq_counter,flip_counter,show_progress=False):

    # Unique list and counts sorted from most to least common
    unique_list = [x[0] for x in vt_counter.most_common()]
    unique_count = [x[1] for x in vt_counter.most_common()]


    # IF READS ARE IN COUNTER:
    if len(unique_list)>0:

        TAGLEN = len(unique_list[0])
        L = len(unique_list)

        ## match_dict stores the first hamming-radius match for each entry (initialized to itself)
        match_dict = np.arange(L)
        ## omit dictionary is a dictionary for each position 0 <= x < TAGLEN
        ## such that when index == index,
        ## omit_dictionary[i] has keys of length TAGLEN -1
        ## derived from the previous entries tag[k] for k < index
        ## by omitting the i-th letter,
        ## and the value is the first index with that key.
        omit_dictionary = [{} for _ in range(TAGLEN)]
        counter = 0

        for index, tag in enumerate(unique_list):
            counter += 1
            if show_progress and counter % 100000 == 0:
                print(counter)
            ## keep a set of values where the i-th omit matches a previous key.
            match       = set()
            for i in range(TAGLEN):
                # print(i)
                newtag = tag[:i] + tag[(i+1):]
                try:
                    ## check if it's already in the dictionary
                    has_match = omit_dictionary[i][newtag]
                    ## if it is, add it to the list of matches
                    match.add(has_match)
                except KeyError:
                    ## if not, there is a error raised
                    ## and we now assign the present index to this entry
                    omit_dictionary[i][newtag] = index
            ## You are assigned to your top match by a sequence of hamming-1 distances
            ## that connect lower counts to higher counts
            matches = list(np.sort(list(match))) + [index]
            match_dict[index]  = match_dict[matches[0]] # select lowest index matching tag - then go look up what value it has stored in the dictionary

        # Now re-process the original counters using match_dict to combine tags
        new_vt_counter     = Counter()
        new_vt_seq_counter = defaultdict(Counter)
        new_flip_counter = defaultdict(double_counter)

        for i in range(len(match_dict)): # go through each unique tag once
            clusterindex = match_dict[i] # get the assignment for this tag
            clustertag = unique_list[clusterindex] # get the tag seq that it was assigned to

            origtag = unique_list[i]
            origcount = unique_count[i] # count of originally unique tag

            new_vt_counter[clustertag] += origcount # update new counter with old count, new tag
            seqs_for_vt = vt_seq_counter[origtag]
            for (r1, r2), readcount in seqs_for_vt.most_common(): # add reads to new counter with new tag
                new_vt_seq_counter[clustertag][(r1,r2)] += readcount

                new_flip_counter[(clustertag,r1,r2)][0] += flip_counter[(origtag,r1,r2)][0]
                new_flip_counter[(clustertag,r1,r2)][1] += flip_counter[(origtag,r1,r2)][1]

        return new_vt_counter, new_vt_seq_counter, unique_list, unique_count, match_dict, new_flip_counter
    else:
        return vt_counter, vt_seq_counter, [] , [], {}, flip_counter





########################################################################
def plot_number_of_reads_per_tag(
                                numreads,
                                numtags,
                                totalreads,
                                totaltags,
                                sample="",
                                filename="number_reads_per_tag.png",
                                logscale=True):

    fig = plt.figure(figsize=(50,10))
    fig.suptitle(sample+"\n"+
                 "Number of Reads Per Tag" + "\n" +
                 "Total Reads: %d" % totalreads + "\n"
                 "Total Tags: %d" % totaltags, fontsize=16)
    width=0.8
    plt.bar(numreads,numtags,0.8,log=logscale)
    #plt.xticks(numreads,numreads)
    plt.xlabel("Reads Per Tag",fontsize=14)
    plt.ylabel("Number of Tags",fontsize=14)
    plt.savefig(filename, dpi=200, facecolor='w', edgecolor='w',
                papertype=None, format=None,
                transparent=False)
    plt.close()


def plot_at_least_x_reads_per_tag(
                                numreads,
                                numtags,
                                totalreads,
                                totaltags,
                                sample="",
                                filename="atleastX_reads_per_tag.png",
                                logscale=True,
                                maxcount=100):
    fig = plt.figure(figsize=(50,10))
    fig.suptitle(sample+"\n"+
                 "Number of Tags with at Least X Reads" + "\n" +
                 "Total Reads: %d" % totalreads + "\n"
                 "Total Tags: %d" % totaltags, fontsize=16)
    x=[]; y=[]
    for xbin in range(1,100):
        x.append(xbin)
        y.append(sum(numtags[numreads>=xbin]))
    plt.bar(x,y,0.8,color='green')
    plt.xlabel("Reads Per Tag",fontsize=14)
    plt.ylabel("Number of Tags with at least X Reads Per Tag",fontsize=14)
    plt.savefig(filename, dpi=200, facecolor='w', edgecolor='w',
                papertype=None, format=None,
                transparent=False)
    plt.close()
########################################################################


def convert_quality_score(qual_ascii):
    qual_numeric = [ord(x)-33 for x in qual_ascii]
    return qual_numeric


########################################################################
def check_tag_structure(vt,structure):
    structureregex=structure
    hits = re.findall('([ACGTN]{3}[AT]{1}[ACGTN]{3}[AT]{1}[ACGTN]{3}[AT]{1}[ACGTN]{3}[AT]{1}[ACGT]{3})',vt)
    return len(hits)>0
    # vt1="ACTTGGTACCGTTTTAAAG" # perfect
    # vt2="ACTGGGTACCGTTTTCAAG" # not
    # structure="NNNWNNNWNNNWNNNWNNN"
