#!/bin/bash

conda config --add channels pkgw
conda config --add channels bcbio
conda config --add channels r
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels defaults

USROPT=""
if [[ -n $1 ]]; then
    USROPT=$1
fi

exists=`conda info --envs|grep -v "#"|awk '{if($1=="j_frisen_1602")print 0}'`
if test "$exists" = "0"; then
    echo "Updating conda environment j_frisen_1602";
    conda install $USROPT -n j_frisen_1602 --file requirements.txt;
else
    echo "Now creating conda environment j_frisen_1602";
    conda create $USROPT -n j_frisen_1602  --file requirements.txt;
fi
echo "To start conda environment type:"
echo "source activate j_frisen_1602"
echo "To stop conda environment type:"
echo "source deactivate [j_frisen_1602]"
# To remove conda environment from system: conda env remove -n
# j_frisen_1602 (convenient if you need to add or update the packages
# installed by conda).
