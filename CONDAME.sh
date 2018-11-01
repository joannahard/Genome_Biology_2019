#!/bin/bash

conda config --add channels pkgw
conda config --add channels bcbio
conda config --add channels r
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels defaults

PROJ="j_frisen_1602_tmp" 

USROPT=""
if [[ -n $1 ]]; then
    USROPT=$1
fi

exists=`conda info --envs|grep -v "#"|awk '{if($1=="$PROJ")print 0}'`
if test "$exists" = "0"; then
    echo "Updating conda environment $PROJ";
    conda install $USROPT -n $PROJ --file requirements.txt;
else
    echo "Now creating conda environment $PROJ";
    conda create $USROPT -n $PROJ  --file requirements.txt;
fi
gatk-register /media/box2/Experiments/Joanna/bin/downloads/GenomeAnalysisTK-3.6.tar.bz2
echo "To start conda environment type:"
echo "source activate $PROJ"
echo "To stop conda environment type:"
echo "source deactivate [$PROJ]"
# To remove conda environment from system: conda env remove -n
# $PROJ (convenient if you need to add or update the packages
# installed by conda).
