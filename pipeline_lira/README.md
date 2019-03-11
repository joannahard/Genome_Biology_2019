# Pipeline for running LiRA.

The LiRA code was downloaded from https://github.com/parklab/LiRA (latest commit c76e9a7cb34af2e44f38bc324f24521ffc1bad59, from Mon Jun 18 16:12:50 2018), but required some changes to work properly with our computational setup and our data; our changes to the code is collected as a patch in the file lira_code_changes/lira.patch. 

The LiRA pipeline requires the latest version of snakemake to handle dynamic rules - use conda "snake2" - requirements file requirements_snake2.txt instead of the conda environment used in the rest of the analysis pipeline. It also needs a conda for all LiRA dependencies with python2 and paths to the lira program needs to be set - requirements file requirements_LiRA.txt. 

Paths were set as:
export PATH=/media/box2/Experiments/Joanna/bin/miniconda3/envs/LiRA/bin:/media/box2/Experiments/Joanna/bin/miniconda3/envs/snake2/bin:$PATH
export PATH=$PATH:/media/box2/Experiments/Joanna/LiRA/LiRA -
export LIRA_DIR=/media/box2/Experiments/Joanna/LiRA/LiRA

There are 2 separate snakemake pipelines, one using dynamic rules for the first steps. Anotherone that takes output from the first one and runs in non-dynamic mode since dynamic was very slow. 

Example scripts for running all steps for all cells in folder scripts.


