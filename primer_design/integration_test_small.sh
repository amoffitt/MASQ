#!/bin/bash

rm -rf design_results
mkdir -p design_results

python select_enzymes_for_snps.py \
    config.primerdesign.example_small.yaml 2>&1 > log.primerdesign_small.txt

diff snp_primer_design_results.primerdesign_example_small.txt \
    design_results/snp_primer_design_results.primerdesign_example_small.*.txt || \
    { echo 'ERROR: design_results/snp_primer_design_results.primerdesign_example_small.*.txt differs from the expected' >&2; exit 1; }

