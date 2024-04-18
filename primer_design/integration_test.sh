#!/bin/bash

rm -rf design_results
mkdir -p design_results

masq_select_enzymes_for_snps \
    config.primerdesign.example.yaml 2>&1 > primerdesign.log

diff snp_primer_design_results.primerdesign_example_updated.txt \
    design_results/snp_primer_design_results.primerdesign_example.*.txt || \
    { echo 'ERROR: design_results/snp_primer_design_results.primerdesign_example.*.txt differs from the expected' >&2; exit 1; }

