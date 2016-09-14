#Introduction

This repository aims to provide the tools to replicate bioinformatics
results obtained in the NBIS long-term support project J_Frisen_1602.
These tools is set to work in a Unix-like operative system (Linux,
Unix, MacOSX, cygwin, etc.) and you will need to work from a
Unix-terminal.

#Results

##Snakemake proceudre for Assembly

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
programs. You should activate the environment before running snakemake

source activate j_frisen_1602

Your terminal prompt will change to show that you are now working in
the j_frisen_1602 conda environment. You can deactivate the
environment by

source deactivate j_frisen_1602

and you will return t your standard terminal environment.

### The snakemake file