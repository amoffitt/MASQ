## Multiplex Accurate Sensitive Quantitation (MASQ) Analysis and Primer Design Design Pipelines


### Set-up
```bash
# Clone workflow into working directory
git clone https://github.com/amoffitt/MASQ /path/to/workdir

# Set-up virtual environment using conda
# Installs required Python and R packages
conda env creat -f environment.yaml
conda activate MASQ

# Edit configuration as needed for MASQ analysis (see below)
vi config.yaml

# Execute MASQ analysis workflow (in dry-run mode)
snakemake -n
```

### MASQ Analysis
The MASQ analysis pipeline is contained in a Snakemake workflow (https://snakemake.readthedocs.io/en/stable/). The workflow is defined by rules in the Snakefile. Individual scripts called from the Snakefile are located in the scripts folder. 

Example input files are included for testing of the workflow installation. Example files include small snippets of FASTQ and BAM files, and a corresponding example locus table. The config.yaml file included here is set-up to run the example files. To execute the workflow on the example, run the snakemake command from the working directory. 

To run on a cluster, a cluster.yaml file can be added to specify parameters specific to your cluster setup, as described in the Snakemake documentation (https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html). 

To run data other than the examples, edit the config.yaml file as necessary to point to FASTQ files, the loci table, WGS bam files, and other sample and run parameters. 


### MASQ Primer Design
The MASQ primer design pipeline is run with the following command: 

```bash
cd primer_design
python select_enzymes_for_snps.py config.primerdesign.example.yaml 2>&1 | tee log.primerdesign.txt
```

Example configuration file and example input SNV lists are included in the primer\_design folder. Enzyme cut site files are included for hg19. Edit the configuration file to point to different SNV files or change the primer design parameters. 

BLAT and Primer3 are installed via the provided conda environment file, but can be installed separately and placed in the path. 
