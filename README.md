# Introduction

This repository aims to provide the tools to  bioinformatics results reported in the manuscript **Conbase: a software for unsupervised discovery of clonal somatic mutations in single cells through read phasing** by Hård et al., submitted to *Genome Biology*.These tools is set to work in a Unix-like operative system (Linux, Unix, MacOSX, cygwin, etc.) from a Unix-terminal.

The main aim is documentation of the assembly and analysis pipelines used. The data used in the manuscript is not included in this repository.  

# Results

## Snakemake procedures for Assembly and simulations

### The conda environment

We use a conda environment to set up the software framework needed to
perform the analyses. You will need to install the miniconda
program, there are good instructions on
http://conda.pydata.org/docs/install/quick.html.

Once, you have miniconda installed, running the script CONDAME.sh:

bash CONDAME.sh

sets up the conda environment 'j_frisen_1602' with the required
programs.
NB! GATK cannot be installed completely from bioconda since it requires
a license. Therefore you have to download GATK (v3.6 or which ever is
specified in requirements.txt) from the website. After activating the
conda environment you need to run gatk-register e.g.:

gatk-register bin/downloads/GenomeAnalysisTK-3.6.tar.bz2

You should activate the environment before running snakemake

source activate j_frisen_1602

Your terminal prompt will change to show that you are now working in
the j_frisen_1602 conda environment. You can deactivate the
environment by

source deactivate j_frisen_1602

and you will return to your standard terminal environment.

### The snakemake files

There are three separate pipelines:

#### The main assembly pipeline
This is defined in `Snakefile` in the main/base directory perform read-mapping, QC and initial variant-calling of single-cell data.

#### The LiRA pipeline
This is defined in `pipeline_lira/Snakefile` and describes how analyses using the program LiRA were run. The folder `pipeline_lira/Snakefile` also contain the relevant. See further the README.md file in the folder `pipeline_lira`.

The LiRA code was downloaded from `https://github.com/parklab/LiRA` (latest commit c76e9a7cb34af2e44f38bc324f24521ffc1bad59, from Mon Jun 18 16:12:50 2018), but required some changes to work properly with our setup and our data (especially our simulation data); our changes to the code is collected as a patch in the file `lira_code_changes/lira_sim.patch`. We plan to share relevant changes (that not specific our simualtions) with the LiRA authors as a pull-request to their code repo.

#### The simulation pipeline
This is defined in `pipeline_simulation/Snakefile` and performs the simulation study from the manuscript, that is, generates synthetic single-cell data based on observed bulk data, runs Conbase, Monovar, LiRA, and SCcaller variant calling on this data and summarizes the results. See further the README.md file in the folder `pipeline_simulation`.

NB! The variant calling programs are not installed automatically and their location needs to be given in the file `config_pipeline.yaml`.

The simulation model description and the scripts used to generate synthetic data can be found in the folder ``simulation``
