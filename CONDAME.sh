#!/bin/bash
echo "Now installing conda environment j_frisen_1602"
conda create -n j_frisen_1602 -c bioconda --file requirements.txt

echo "To start conda environment type:"
echo "source activate j_frisen_1602"
echo "To stop conda environment type:"
echo "source deactivate [j_frisen_1602]"
# To remove conda environment from system: conda env remove -n
# j_frisen_1602 (convenient if you need to add or update the packages
# installed by conda).
