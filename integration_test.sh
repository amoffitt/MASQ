#!/bin/bash


./prepare_example_references.sh

rm -rf example_outputs
rm -rf test_example1
rm -rf test_example2
rm -rf combined
rm -rf .snakemake

tar zxf example_outputs.tar.gz

snakemake -j 10

diff -r -x logs -x fastqc \
    example_outputs/test_example1 \
    test_example1 || { echo 'test_example1 differs from the expected' >&2; exit 1; }
diff -r  -x logs -x fastqc \
    example_outputs/test_example2 \
    test_example2 || { echo 'test_example2 differs from the expected' >&2; exit 1; }

