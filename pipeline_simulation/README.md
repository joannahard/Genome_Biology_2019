# Pipeline for simulation

Performs the simulation study in the manuscript, that is, generates synthetic single-cell data based on observed bulk data, runs Conbase, Monovar, LiRA, and SCcaller variant calling on the data and summarizes the results. 

NB! The variant calling programs are not installed automatically and their location needs to be given in the file config_pipeline.yaml.

The simulation model description and the scripts used to generate synthetic data can be found in the folder simulation




To run the pipeline, the conda "snake2" has all requirements - requirements for conda env in requirements_snake2.txt. In addition, conda environment for LiRA is required, more info in folder pipeline_lira. Also, for running SCcaller, another conda env was used, see requirements_sccaller.txt. 

# run full pipeline to merged stats matrix:
source activate snake2
snakemake -p -j 80 merge_stats_all

# the 2 configfiles are:
- config_pipeline.yaml - simulated depth
- config_pipeline_flat.yaml - flat depth at 30

# some target rules:
- all - will exectute all the steps for conbase,sccaller and monovar. Lira is not included in the all rule.
- all_conbase - run all steps of conbase
- all_gatk - preprocessing with gatk for lira
- all_sim - create all simulation datasets.
- all_lira - run all steps of lira
- all_monovar
- all_sccaller


NB! summary plotting script is not included in pipeline, they are included in folder scripts/R/ and were run manually.