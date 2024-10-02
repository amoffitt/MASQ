#!/usr/bin/env bash

pip install -e . > pipeline.log 2>&1

./prepare_example_references.sh

./integration_cleanup.sh

mkdir -p example_outputs
tar zxf  masq_full_example_outputs_sample1.tar.gz -C example_outputs
tar zxf  masq_full_example_outputs_sample2.tar.gz -C example_outputs

cd examples/masq_example
snakemake -j > pipeline.log 2>&1
cd -

diff -r -x logs -x fastqc \
    example_outputs/sample1_cellfree \
    examples/masq_example/sample1_cellfree \
    || { echo 'sample1_cellfree differs from the expected' >&2; exit 1; }
diff -r  -x logs -x fastqc \
    example_outputs/sample2_cells \
    examples/masq_example/sample2_cells \
    || { echo 'sample2_cells differs from the expected' >&2; exit 1; }

