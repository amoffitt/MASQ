#!/bin/bash


mkdir -p references; cd references
mkdir -p hg19; cd hg19

if [[ ! -f hg19.fa ]]; then
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
    gunzip hg19.fa.gz
fi

if [[ ! -f hg19.fa.fai ]]; then
    samtools faidx hg19.fa
fi

cd ../..
