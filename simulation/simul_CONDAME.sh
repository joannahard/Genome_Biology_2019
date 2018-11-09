#!/bin/bash

conda config --add channels clinicalgraphics
#conda config --add channels bcbio
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

USROPT=""
if [[ -n $1 ]]; then
    USROPT=$1
fi

PROJ="simul_j_frisen_1602"
REQ="simul_requirements.txt"

exists=`conda info --envs|grep -v "#"|awk -v var="$PROJ" '{if($1==var)print 0}'`
if test "$exists" = "0"; then
    echo "Updating conda environment $PROJ";
    conda install $USROPT -n $PROJ --file $REQ;
else
    echo "Now creating conda environment $PROJ";
    conda create $USROPT -n $PROJ  --file $REQ;
fi

echo "To start conda environment type:"
echo "source activate $PROJ"
echo "To stop conda environment type:"
echo "source deactivate [$PROJ]"
# To remove conda environment from system: conda env remove -n
# $PROJ (convenient if you need to add or update the packages
# installed by conda).
