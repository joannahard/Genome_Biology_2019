import pandas as pd
import os, sys, math


from argparse import ArgumentParser
parser = ArgumentParser(description='Calculate qc stats based on a gene expression matrix')

# [Required input]
parser.add_argument('-m', '--monovar_file', metavar='monovar_file', help="Monovar vcf file", required=True)
parser.add_argument('-c', '--conbase_file', metavar='conbase_file', help="Conbase tsv stats file", required=True)
parser.add_argument('-s', '--sim_file', metavar='sim_file', help="Simulation stats json file", required=True)
parser.add_argument('-o', '--outfile', metavar='outfile', help='Output file', required=True)
# [Optional input]
parser.add_argument('--monovar_dp',metavar='monovar_dp', help="DP cutoff for parsing monovar", type =int, default = 5)
parser.add_argument('--monovar_gq',metavar='monovar_gq', help="GQ cutoff for parsing monovar", type =int, default = 0)
parser.add_argument('--monovar_nmut',metavar='monovar_nmut', help="Required number of cells with called mutation in monovar", type =int, default = 2)
args = parser.parse_args()

#################################################3
# check inputs:

for file in [args.monovar_file, args.conbase_file, args.sim_file]:
    if not os.path.exists(file):
        print("Error! No such file "+ file)
        sys.exit(1)


monovar_dp = args.monovar_dp
monovar_gq = args.monovar_gq
monovar_nmut = args.monovar_nmut


#################################################3
# some functions:

def check_mut(sdict):
    r = sdict['Ref'].split('_')
    m = sdict['Mut'].split('_')
    if r[1] == m[1]:
        return(False)
    else:
        return(True)

def parse_monovar(stats):
    s = stats.split(':') # will have entries: ['GT','AD','DP','GQ','PL']
    gt = s[0]
    if len(s) == 1 and gt == "./.":
        return(-1)
    if int(s[2]) < monovar_dp:
        return(-1)
    if int(s[3]) < monovar_gq:
        return(-1)
    if gt == "0/0":
        return(0)
    else:
        return(1)

#################################################3
# read as pandas df
sim_data = pd.read_json(args.sim_file)
sim_data = sim_data.T # transpose

# extend the ado data and states data
ado = sim_data.ADO.apply(pd.Series)
states = sim_data.states.apply(pd.Series)

# for each cell the state tells if the cell is mutated or not.
# e.g.: cell0_state     {'Ref': 'C_T', 'Mut': 'T_G'}    True - has mut in clone1
# cell10_state    {'Ref': 'C_G', 'Mut': 'T_G'}   False - has no mut in clone2
# need to parse the dict and get assignment

has_mut = states.apply(lambda x: x.apply(check_mut), axis=1)
has_mut.columns = [str(col) + "_mut" for col in states.columns]
sim_data = pd.concat([sim_data.drop(["ADO","states"], axis = 1), ado, has_mut], axis=1)

######################################################################    
# read conbase data
conbase_raw = pd.read_csv(args.conbase_file, sep= "\t")

# reads in extra column "Unnamed: 25"
conbase_raw = conbase_raw[conbase_raw.columns.drop(list(conbase_raw.filter(regex='Unnamed')))]    

# also prints rows with CONFLICT - remove all such rows.
has_conflict = conbase_raw.apply(lambda x: any(x[5:].str.contains("CONFLICT")), axis=1)
if conbase_raw.shape[0] > 0 and any(has_conflict):
    conbase_raw = conbase_raw.drop(conbase_raw.index[has_conflict])

cb_sites = conbase_raw.CHROM.map(str) + ":" + conbase_raw.POS.map(str)
conbase_raw.index = cb_sites
#filter for sites in sim_data
conbase_raw = conbase_raw[conbase_raw.index.isin(sim_data.index)]


cb_cells = [x.replace(":DP","") for x in list(conbase_raw.columns)[5:]]
cb_cells = [x.replace("sim_","") for x in cb_cells]


# convert to 1,0 or -1 (if non-informative)
conbase = conbase_raw.iloc[:,5:]
conbase.columns = cb_cells
conbase.replace({r'HOMO-C.:.+': '0'}, regex=True, inplace = True)
conbase.replace({r'HET-C.:.+': '1'}, regex=True, inplace = True)
conbase.replace({r'NOT-INFORMATIVE:.+': '-1' }, regex=True, inplace = True)
conbase.replace({r'ZERO-READS:0': '-1' }, regex=True, inplace = True)


# convert to int or na
conbase = conbase.astype(int)


######################################################################
# read monovar data

monovar_raw = pd.read_table(args.monovar_file,sep="\t",header=19,index_col=False)
monovar_raw.index = monovar_raw['#CHROM'].map(str) + ":" + monovar_raw['POS'].map(str).tolist()
# filter for sites in sim_data
monovar_raw = monovar_raw[monovar_raw.index.isin(sim_data.index)]

# set all sites with DP < monovar_dp or GQ < monovar_gq  to NA
# set all with 0/0 to 0
# set all with 1/0 or 1/1 to 1
monovar = monovar_raw[cb_cells].apply(lambda x: x.apply(parse_monovar), axis=1)

# only keep cells with >= nmut cells with 1
if monovar.shape[0] > 0:
    monovar_counts = monovar.apply(pd.Series.value_counts, axis=1)
    monovar = monovar[monovar_counts[1] >= monovar_nmut]


# OBS! almost all sites with 0/1 get filtered due to GQ < 5

######################################################################
# get resuls per site/cell

# consider results as:
# T-het - has correctly predicted genotype ALT when has Mut.
# T-hom - has correctly predicted genotype REF when no Mut
# F-het - has wrong genotype - has Mut but predicts REF
# F-hom - has wrong genotype - has no Mut but predicts ALT
# NA - has not predicted anything.

# also add on class - wEAL

def summarize(pred,sim_data):
    data = pd.DataFrame(index=sim_data.index, columns = ['EAL', 'wADO', 'noADO', 'wMut','noMut','T-het','T-hom','F-het','F-hom','NA','T-het_wEAL','T-hom_wEAL','F-het_wEAL','F-hom_wEAL','NA_wEAL'])
    data[data.columns[5:]] = 0
    data["EAL"] = sim_data["EAL"]
    mut_columns = [x + '_mut' for x in cb_cells]

    for site in data.index:
        vc = sim_data.ix[site][cb_cells].value_counts()
        data.ix[site,'wADO'] = vc[True] if True in vc else 0
        data.ix[site,'noADO'] = vc[False] if False in vc else 0
    
        vc = sim_data.ix[site][mut_columns].value_counts()
        data.ix[site,'wMut']= vc[True] if True in vc else 0
        data.ix[site,'noMut']= vc[False] if False in vc else 0

        for cell in cb_cells:
            mut = sim_data.ix[site][cell + "_mut"]  # mut=True has mut, mut=False no mut
            if site in pred.index:
                result = ""
                cb = pred.ix[site][cell]
                if cb == -1:
                    result = "NA"
                elif mut and cb == 1: # has mut and cb calls 1
                    result = "T-het"
                elif mut and cb == 0:
                    result = "F-het"
                elif (not mut) and cb == 0:
                    result = "T-hom"
                elif (not mut) and cb == 1:
                    result = "F-het"
                else:
                    print("Error! inconsistent predictions for " + site + " : " + cell)
                    print("pred_data")
                    print(pred.ix[site])
                    print("sim_data")
                    print(sim_data.ix[site])
                    sys.exit(1)
            else:
                result = "NA"
            if sim_data.ix[site]["EAL"]:
                result = result + "_wEAL"
            data.ix[site,result] += 1
    #counts = data[data.columns[5:]].apply(pd.Series.value_counts, axis=0)
    #counts.fillna(0, inplace = True)
    #counts = counts.astype(int)
    #counts.sort_index(axis=1, inplace=True)
    counts = data[data.columns[5:]].sum()
    return(data,counts)

cb_data,cb_counts = summarize(conbase, sim_data)
#print("conbase")
#print(cb_counts)

mv_data, mv_counts = summarize(monovar, sim_data)
#print("monovar")
#print(mv_counts)

out = pd.DataFrame()
out["conbase"]=cb_counts
out["monovar"]=mv_counts
#print(out)
out.to_csv(args.outfile)
