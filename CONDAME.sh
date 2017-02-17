#!/bin/bash


exists=`conda info --envs|grep -v "#"|awk '{if($1=="j_frisen_1602")print 0}'`
if test "$exists" = "0"; then
    echo "Updating conda environemnt j_frisen_1602";
    conda install -n j_frisen_1602 -c bioconda -c pkgw --file requirements.txt;
else
    echo "Now creating conda environment j_frisen_1602";
    conda create -n j_frisen_1602 -c bioconda -c pkgw --file requirements.txt;
fi
echo "To start conda environment type:"
echo "source activate j_frisen_1602"
echo "To stop conda environment type:"
echo "source deactivate [j_frisen_1602]"
# To remove conda environment from system: conda env remove -n
# j_frisen_1602 (convenient if you need to add or update the packages
# installed by conda).
