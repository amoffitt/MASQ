#!/usr/bin/env bash

pip install -e . > pipeline.log 2>&1

./prepare_example_references.sh

./integration_cleanup.sh

mkdir -p example_outputs
tar zxf masq_example_outputs_noqc_double_counter_qcfiltered.tar.gz -C example_outputs

ln -s config.masq.yaml config.yaml

snakemake -j 5 > pipeline.log 2>&1

diff -r -x logs -x fastqc \
    example_outputs/test_example1 \
    test_example1 || { echo 'test_example1 differs from the expected' >&2; exit 1; }
diff -r  -x logs -x fastqc \
    example_outputs/test_example2 \
    test_example2 || { echo 'test_example2 differs from the expected' >&2; exit 1; }

