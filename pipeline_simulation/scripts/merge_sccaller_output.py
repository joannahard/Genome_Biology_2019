import pandas as pd
import numpy as np
import os, sys, math

# Will take ouput from sccaller for all cells with one simulation setting and parse into one matrix.
# for each site put 0/1 if PASS for sccaller.

from argparse import ArgumentParser
parser = ArgumentParser(description='Calculate qc stats based on a gene expression matrix')

# [Required input]
parser.add_argument('-i', '--infiles', metavar='infiles', nargs = "+", help="SCCaller files", required=True)
parser.add_argument('-o', '--outfile', metavar='outfile', help='Output file', required=True)
args = parser.parse_args()


#################################################3
# read as pandas df

#args.infiles = ["data/sim_snv1_eal0.5_ado0.5/sccaller/output/sim_cell0.sincell_filtered0.01.bed","data/sim_snv1_eal0.5_ado0.5/sccaller/output/sim_cell10.sincell_filtered0.01.bed","data/sim_snv1_eal0.5_ado0.5/sccaller/output/sim_cell11.sincell_filtered0.01.bed","data/sim_snv1_eal0.5_ado0.5/sccaller/output/sim_cell12.sincell_filtered0.01.bed","data/sim_snv1_eal0.5_ado0.5/sccaller/output/sim_cell13.sincell_filtered0.01.bed","data/sim_snv1_eal0.5_ado0.5/sccaller/output/sim_cell14.sincell_filtered0.01.bed","data/sim_snv1_eal0.5_ado0.5/sccaller/output/sim_cell15.sincell_filtered0.01.bed","data/sim_snv1_eal0.5_ado0.5/sccaller/output/sim_cell16.sincell_filtered0.01.bed","data/sim_snv1_eal0.5_ado0.5/sccaller/output/sim_cell17.sincell_filtered0.01.bed","data/sim_snv1_eal0.5_ado0.5/sccaller/output/sim_cell18.sincell_filtered0.01.bed","data/sim_snv1_eal0.5_ado0.5/sccaller/output/sim_cell19.sincell_filtered0.01.bed","data/sim_snv1_eal0.5_ado0.5/sccaller/output/sim_cell1.sincell_filtered0.01.bed","data/sim_snv1_eal0.5_ado0.5/sccaller/output/sim_cell2.sincell_filtered0.01.bed","data/sim_snv1_eal0.5_ado0.5/sccaller/output/sim_cell3.sincell_filtered0.01.bed","data/sim_snv1_eal0.5_ado0.5/sccaller/output/sim_cell4.sincell_filtered0.01.bed","data/sim_snv1_eal0.5_ado0.5/sccaller/output/sim_cell5.sincell_filtered0.01.bed","data/sim_snv1_eal0.5_ado0.5/sccaller/output/sim_cell6.sincell_filtered0.01.bed","data/sim_snv1_eal0.5_ado0.5/sccaller/output/sim_cell7.sincell_filtered0.01.bed","data/sim_snv1_eal0.5_ado0.5/sccaller/output/sim_cell8.sincell_filtered0.01.bed","data/sim_snv1_eal0.5_ado0.5/sccaller/output/sim_cell9.sincell_filtered0.01.bed"]
#args.outfile = "temp/scc_stats.csv"


all = pd.DataFrame()

for file in args.infiles:
    cell = os.path.split(file)[1].replace(".","_").split("_")[1]
    file_size = os.path.getsize(file)
    if file_size > 0:
        d = pd.read_csv(file, sep="\t", header=None)
        sites = d[0].map(str) + ":" + d[1].map(str)
        d2 = pd.DataFrame(index=sites)
        d2[cell] = 1
        all = pd.concat([all,d2], axis=1)
    else: # if empty file, add a column with NaN
        all[cell] = np.nan
        
all.to_csv(args.outfile, sep=",")
