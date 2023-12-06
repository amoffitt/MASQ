'''
Helper script to convert genome FASTA into python dictionary and store as pickle file
'''

import os
import sys
import pickle
import sys

def make_ref_genome_pickle(ref_fasta, ref_pickle):

    if not os.path.exists(ref_pickle):
        print("Pickle file (%s) does not exist. Creating from FASTA (%s)." % (ref_pickle, ref_fasta))
        seqdic = dict()
        chr = ""
        with open(ref_fasta, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    print(line)
                    # store previous chromosome
                    if len(chr) > 0:
                        seqdic[chr] = "".join(seqs)
                    # start new chromosome
                    chr = line.strip()[1:]
                    seqs = []
                else:
                    seqs.append(line.strip().upper())
            # store last chromosome
            seqdic[chr] = "".join(seqs)

        pickle.dump(seqdic, open(ref_pickle, "wb"))


if __name__=="__main__":
    ref_fasta=sys.argv[1]
    ref_pickle=sys.argv[2]
    make_ref_genome_pickle(ref_fasta, ref_pickle)
