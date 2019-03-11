# Introduction

This repository aims to provide the tools to  bioinformatics results reported in the manuscript **Conbase: a software for unsupervised discovery of clonal somatic mutations in single cells through read phasing** by Hård et al. (*Genome Biology*). These tools is set to work in a Unix-like operative system (Linux, Unix, MacOSX, cygwin, etc.) from a Unix-terminal.

The main aim is documentation of the data analysis pipelines used. The data used in the manuscript is available on https://ega-archive.org/studies/EGAS00001003108.  

# Results

## Snakemake procedures for data processing and simulations

### The conda environments

We have used conda environments to set up the software framework needed to
perform the analyses (https://docs.conda.io). The required conda-packages are listed in pipeline-spcific (see below) `requirements.yaml/txt` files  

NB! GATK cannot be installed completely from bioconda since it requires
a license. Therefore you have to download GATK (v3.6 or which ever is
specified in requirements.txt) from the website. After activating the
conda environment you need to run `gatk-register` and refer to the downloaded GATK.

### The snakemake files

There are three separate pipelines:

#### The main variant-calling pipeline
This is defined in `Snakefile` in the main/base directory perform read-mapping, QC and initial variant-calling of single-cell data. Conda requirements can be found in `requirements.yaml`.

This pipeline was used for the analyses of both experimental single-cell datasets (human T-cell and fibroblasts, respectively) in the manuscript. However, different configuration files were used for the analyses: `config_Tcell.yaml`, `config_Fibs.yaml`, respectively.

#### The LiRA pipeline
This is defined in `pipeline_lira/Snakefile` and describes how analyses using the program LiRA were run. The folder `pipeline_lira/Snakefile` also functions as the relevant documentation. See further the README.md file in the folder `pipeline_lira`.

#### The simulation pipeline
This is defined in `pipeline_simulation/Snakefile` and performs the simulation study from the manuscript, that is, generates synthetic single-cell data based on observed bulk data, runs Conbase, Monovar, LiRA, and SCcaller variant calling on this data and summarizes the results. See further the README.md file in the folder `pipeline_simulation`.

The simulation model description and the scripts used to generate synthetic data can be found in the folder ``simulation``
