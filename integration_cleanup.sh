#!/usr/bin/env bash

rm -rf example_outputs
rm -rf test_example1
rm -rf test_example2
rm -rf sample1_blood
rm -rf sample2_tumor
rm -rf sample1_cellfree
rm -rf sample2_cells
rm -rf combined
rm -rf .snakemake
rm -rf config.yaml

rm -rf examples/masq_example/.snakemake
rm -rf examples/masq_example/sample1_cellfree
rm -rf examples/masq_example/sample2_cells
rm -rf examples/masq_example/combined
rm -rf examples/masq_example/pipeline.log

rm -rf examples/pcr_example/.snakemake
rm -rf examples/pcr_example/sample1_blood
rm -rf examples/pcr_example/sample2_tumor
rm -rf examples/pcr_example/combined
rm -rf examples/pcr_example/pipeline.log
