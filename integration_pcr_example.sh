#!/usr/bin/env bash

pip install -e . > pipeline.log 2>&1

./prepare_example_references.sh

./integration_cleanup.sh

mkdir -p examples_output/outputs
tar zxf examples_output/pcr_example_outputs_fixed.tar.gz -C examples_output/outputs


cd examples/pcr_example
snakemake -j > pipeline.log 2>&1
cd -

diff -r -x logs -x fastqc \
    examples/pcr_example/sample1_blood \
    examples_output/outputs/sample1_blood \
    || { echo 'sample1_blood differs from the expected' >&2; exit 1; }
diff -r  -x logs -x fastqc \
    examples/pcr_example/sample2_tumor \
    examples_output/outputs/sample2_tumor \
    || { echo 'sample2_tumor differs from the expected' >&2; exit 1; }

