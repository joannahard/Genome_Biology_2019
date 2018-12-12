import pandas as pd
import os, sys, math


from argparse import ArgumentParser
parser = ArgumentParser(description='Merge all stats from lira,conbase,monovar & sccaller')

# [Required input]
parser.add_argument('-o', '--outprefix', metavar='outprefix', help='Outprefix', required=True)
parser.add_argument('-i', '--infiles', metavar='infiles', help='Input stats files', required=True, nargs="+")

args = parser.parse_args()

#################################################3
# check inputs:

outCB = pd.DataFrame()
outMV = pd.DataFrame()
outSCC = pd.DataFrame()
outL = pd.DataFrame()

for infile in args.infiles:
    if not os.path.exists(infile):
        print("Error! No such file "+ infile)
        sys.exit(1)
    data = pd.read_csv(infile, sep = ",", index_col = 0)
    sim_name = infile.split("/")[1]
    settings = sim_name.split("_")[1:]  # ['snv1', 'eal0.5', 'ado0.5']
    snv = float(settings[0].replace("snv",""))
    eal = float(settings[1].replace("eal",""))
    ado = float(settings[2].replace("ado",""))        
    data.loc["fSNV"] = [snv,snv,snv,snv]
    data.loc["fEAL"] = [eal,eal,eal,eal]        
    data.loc["fADO"] = [ado,ado,ado,ado]
    outCB[sim_name] = data["conbase"]
    outMV[sim_name] = data["monovar"]
    outL[sim_name] = data["lira"]
    outSCC[sim_name] = data["sccaller"]        

#print(out)
outCB.to_csv(args.outprefix + "_conbase.csv")
print("Conbase stats written to "+ args.outprefix + "_conbase.csv")

outMV.to_csv(args.outprefix + "_monovar.csv")
print("Monovar stats written to "+ args.outprefix + "_monovar.csv")

outL.to_csv(args.outprefix + "_lira.csv")
print("Lira stats written to "+ args.outprefix + "_lira.csv")

outSCC.to_csv(args.outprefix + "_sccaller.csv")
print("SCCaller stats written to "+ args.outprefix + "_sccaller.csv")


