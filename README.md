# Introduction

This repository aims to provide the tools to replicate bioinformatics
results obtained in the NBIS long-term support project J_Frisen_1602.
These tools is set to work in a Unix-like operative system (Linux,
Unix, MacOSX, cygwin, etc.) and you will need to work from a
Unix-terminal.

# Results

## Snakemake proceudre for Assembly

### The conda environment

We use a conda environment to set up the software framework needed to
perform the analyses. You will need to install the miniconda
program. On UPPMAX, miniconda is available as a module :

module add miniconda3

For other systems, there are good instructions on
http://conda.pydata.org/docs/install/quick.html.

Once, you have miniconda installed, running the script CONDAME.sh:

bash CONDAME.sh

sets up the conda environment 'j_frisen_1602' with the required
programs.
NB! GATK cannot be installed completely from bioconda since it requires
a license. Therefore you have to download GATK (v3.6 or which ever is
specified in requirements.txt) from the website. After activating the
conda environment you need to run gatk-register e.g.:

gatk-register /media/box2/Experiments/Joanna/bin/downloads/GenomeAnalysisTK-3.6.tar.bz2

You should activate the environment before running snakemake

source activate j_frisen_1602

Your terminal prompt will change to show that you are now working in
the j_frisen_1602 conda environment. You can deactivate the
environment by

source deactivate j_frisen_1602

and you will return to your standard terminal environment.

### The snakemake files

There are three separate pipelines:

- The main pipeline, defined in `Snakefile` in the main/base directory perform read-mapping, QC and initial variant-calling of sc data.
- The second pipeline, defined in `pipeline_lira/Snakefile` describes how LiRA analyses were run
- The last pipeline, defined in `pipeline_simulation/Snakefile` performs a simulation study, generating synthetic sc data based on observed bulk data, runs Conbase, Monovar, LiRA, and SCcaller variant calling on this data and summarizes the results. NB! The variant calling programs are not installed automatically and their location needs to be given in the file `config_pipeline.yaml`.
