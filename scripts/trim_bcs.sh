#!/bin/bash

# In-line barcode trimming script

threads=$1
outdir=$2
fastq1=$3
fastq2=$4
sample_bc_list=$5

barcode_program="scripts/trimBarcodeFragments"

$barcode_program -nThreads=$threads -dirPrefix=$outdir -inputFiles=$fastq1,$fastq2 $sample_bc_list

