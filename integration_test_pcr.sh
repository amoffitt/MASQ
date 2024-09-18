#!/usr/bin/env bash

pip install -e . > pipeline.log 2>&1

./prepare_example_references.sh

./integration_cleanup.sh

mkdir -p example_outputs
tar zxf pcr_example_outputs_double_counter_qcfiltered.tar.gz -C example_outputs

ln -s config.standardPCR.yaml config.yaml

snakemake -j 5 > pipeline.log 2>&1

diff -r -x logs -x fastqc \
    example_outputs/sample1_blood \
    sample1_blood || { echo 'sample1_blood differs from the expected' >&2; exit 1; }
diff -r  -x logs -x fastqc \
    example_outputs/sample2_tumor \
    sample2_tumor || { echo 'sample2_tumor differs from the expected' >&2; exit 1; }

