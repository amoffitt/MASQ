## Multiplex Accurate Sensitive Quantitation (MASQ) Analysis and Primer Design Design Pipelines


### Set-up software environment

Clone the `MASQ` repository

```bash
git clone https://github.com/amoffitt/MASQ /path/to/workdir
```

Create a MASQ conda environment using

```bash
conda env create -f environment.yaml
```

and activate it

```bash
conda activate MASQ
```

In your `MASQ` conda environment install the `masq` package. To this end from
the MASQ working directory run:

```bash
pip install -e .
```

### Setup the reference genome

To work with `MASQ` package you will need a reference genome.

To download and prepare the reference genome used by MASQ package and scripts
run the following command:

```bash
bash ./prepare_example_references.sh
```
This command creates a subdirectory `references/hg19` that contains
`hg19` reference genome and a fasta index:

```
references/
└── hg19
    ├── hg19.fa
    └── hg19.fa.fai
```

Alternatively you can download a reference genome and run
`samtools faidx` over the genome fasta file to create a fasta index file.

### MASQ Analysis

The MASQ analysis pipeline is contained in a Snakemake workflow 
(https://snakemake.readthedocs.io/en/stable/). The workflow is defined by 
rules in the Snakefile. Individual scripts called from the Snakefile are 
located in the scripts folder. 

The source code includes two examples.
Example input files are included for testing of the workflow installation. 
Example files include small snippets of FASTQ and BAM files, and a 
corresponding example locus table. The configuration files are included 
to run the example analysis. 

### Setup MASQ project and run the MASQ pipeline

To setup a MASQ project you should choose a directory and run

```
masq_project
```

command. This command will create a `Snakefile` and `config.yaml` file inside
the directory you are running it. 

After that you need to edit the `config.yaml` configuration file to describe
your project. You can check examples configuration files (see `examples` 
subdirectory) for examples how to configure your project.

When the project is configured you can run the MASQ pipeline using:

```
snakemake -j
```


#### Run `MASQ` example

You can look into `examples/masq_example` for an example how to configure MASQ 
pipeline. 

If you enter into `examples/masq_example` directory you can run the MASQ
piplene using:

```bash
snakemake -j
```

#### Run `PCR` example

The `PCR` example is setup in `examples/pcr_example` subdirectory.

If you enter into `examples/pcr_example` directory you can run the MASQ pipeline
using:

```bash
snakemake -j
```

### Run MASQ on cluster

To run on a cluster, a cluster.yaml file can be added to specify parameters specific to your cluster setup, as described in the Snakemake documentation (https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html). 

To run data other than the examples, edit the config.yaml file as necessary to point to FASTQ files, the loci table, WGS bam files, matching reference files, and other sample and run parameters. 


### MASQ Primer Design
The MASQ primer design script is run with the command `masq_select_enzymes_for_snps`.
Using the `-h` or `--help` option the help message is shown:

```bash
masq_select_enzymes_for_snps -h
```

```
usage: masq_select_enzymes_for_snps [-h] config_file

Selects enzymes for SNPs and designs primers

positional arguments:
  config_file  primer design config file

options:
  -h, --help   show this help message and exit
```

### MASQ Primer Design examples

Example configuration files and example input SNV lists are included in the
`primer_design` folder. Enzyme cut site files are included for hg19. Edit the 
configuration file to point to different SNV files or change the primer 
design parameters. 

The source code includes two primer design examples. These examples are
located in the `primer_design` subdirectory of the MASQ source tree. 

Theses
examples are defined by the following configuration files:

```
primer_design/config.primerdesign.example_small.yaml
primer_design/config.primerdesign.example.yaml
```

To run the *small* example use:

```bash
cd primer_design
masq_select_enzymes_for_snps config.primerdesign.example_small.yaml
```

To run the other example use:

```bash
masq_select_enzymes_for_snps config.primerdesign.example.yaml
```

### MASQ developement

Update the MASQ environment with dependencies needed for development:

```bash
mamba env update -f dev-environment.yaml
```

In the activated `MASQ` development environment and from the `MASQ` root
directory run:

```bash
pip install -e .
```

to install an editable version of the `masq` package.
