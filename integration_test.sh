#!/bin/bash

pip install -e .

./prepare_example_references.sh

rm -rf example_outputs
rm -rf test_example1
rm -rf test_example2
rm -rf sample1_blood
rm -rf sample2_tumor
rm -rf combined
rm -rf .snakemake
rm -rf config.yaml

tar zxf masq_example_outputs.tar.gz

ln -s config.masq.yaml config.yaml

snakemake -j 5

diff -r -x logs -x fastqc \
    example_outputs/test_example1 \
    test_example1 || { echo 'test_example1 differs from the expected' >&2; exit 1; }
diff -r  -x logs -x fastqc \
    example_outputs/test_example2 \
    test_example2 || { echo 'test_example2 differs from the expected' >&2; exit 1; }

