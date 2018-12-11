import pandas as pd
import os, sys, math

# Will take ouput from lira for all cells with one simulation setting and parse into one matrix.
# for each site put 0/1 if PASS 

from argparse import ArgumentParser
parser = ArgumentParser(description='Merge lira output')

# [Required input]
parser.add_argument('-i', '--infiles', metavar='infiles', nargs = "+", help="SCCaller files", required=True)
parser.add_argument('-o', '--outfile', metavar='outfile', help='Output file', required=True)
args = parser.parse_args()

#################################################3
# read as pandas df

#infiles = "data/sim_snv1_eal0_ado0.1/lira/lira_output/cell0/varcall_bulk/out.vcf.gz","data/sim_snv1_eal0_ado0.1/lira/lira_output/cell10/varcall_bulk/out.vcf.gz","data/sim_snv1_eal0_ado0.1/lira/lira_output/cell11/varcall_bulk/out.vcf.gz","data/sim_snv1_eal0_ado0.1/lira/lira_output/cell12/varcall_bulk/out.vcf.gz","data/sim_snv1_eal0_ado0.1/lira/lira_output/cell13/varcall_bulk/out.vcf.gz","data/sim_snv1_eal0_ado0.1/lira/lira_output/cell14/varcall_bulk/out.vcf.gz","data/sim_snv1_eal0_ado0.1/lira/lira_output/cell15/varcall_bulk/out.vcf.gz","data/sim_snv1_eal0_ado0.1/lira/lira_output/cell16/varcall_bulk/out.vcf.gz","data/sim_snv1_eal0_ado0.1/lira/lira_output/cell17/varcall_bulk/out.vcf.gz","data/sim_snv1_eal0_ado0.1/lira/lira_output/cell18/varcall_bulk/out.vcf.gz","data/sim_snv1_eal0_ado0.1/lira/lira_output/cell19/varcall_bulk/out.vcf.gz","data/sim_snv1_eal0_ado0.1/lira/lira_output/cell1/varcall_bulk/out.vcf.gz","data/sim_snv1_eal0_ado0.1/lira/lira_output/cell2/varcall_bulk/out.vcf.gz","data/sim_snv1_eal0_ado0.1/lira/lira_output/cell3/varcall_bulk/out.vcf.gz","data/sim_snv1_eal0_ado0.1/lira/lira_output/cell4/varcall_bulk/out.vcf.gz","data/sim_snv1_eal0_ado0.1/lira/lira_output/cell5/varcall_bulk/out.vcf.gz","data/sim_snv1_eal0_ado0.1/lira/lira_output/cell6/varcall_bulk/out.vcf.gz","data/sim_snv1_eal0_ado0.1/lira/lira_output/cell7/varcall_bulk/out.vcf.gz","data/sim_snv1_eal0_ado0.1/lira/lira_output/cell8/varcall_bulk/out.vcf.gz","data/sim_snv1_eal0_ado0.1/lira/lira_output/cell9/varcall_bulk/out.vcf.gz"

#outfile = "temp/lira_all.csv"

all = pd.DataFrame()

for file in args.infiles:
    cell = file.split("/")[4]
    d = pd.read_csv(file, sep="\t", header=None, comment="#", index_col=None)
    # only keep rows with PASS in column 6
    d = d[d[6].str.contains("PASS")]
    sites = d[0].map(str) + ":" + d[1].map(str)
    d2 = pd.DataFrame(index=sites)
    d2[cell] = 1
    all = pd.concat([all,d2], axis=1)

all.to_csv(args.outfile, sep=",")
