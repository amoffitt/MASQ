#!/bin/bash

mkdir references; cd references
mkdir hg19; cd hg19
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
python ../../scripts/make_ref_genome_pickle.py hg19.fa hg19.seq_dic.cpickle
cd ../..
