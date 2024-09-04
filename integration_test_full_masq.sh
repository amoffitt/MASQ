#!/usr/bin/env bash

pip install -e . > pipeline.log 2>&1

./prepare_example_references.sh

./integration_cleanup.sh

mkdir -p example_outputs
tar zxf masq_full_example_outputs_double_counter.tar.gz -C example_outputs

ln -s config.masq.full_example.yaml config.yaml

snakemake -j > pipeline.log 2>&1

diff -r -x logs -x fastqc \
    example_outputs/sample1_cellfree \
    sample1_cellfree || { echo 'sample1_cellfree differs from the expected' >&2; exit 1; }
diff -r  -x logs -x fastqc \
    example_outputs/sample2_cells \
    sample2_cells || { echo 'sample2_cells differs from the expected' >&2; exit 1; }

