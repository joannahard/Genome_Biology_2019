#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys, os, re

def getMapperInput(config, wc, read):
    if len(config["metafiles"][wc.sample]) == 1:
        name = config["metafiles"][wc.sample][0][0:-2]
        # This is the case when the input file has no extension (single experiemtn)
        # Right now, we allow the user to decide what extension shoul be added to
        # the output file. If we want to force the extension to be '1', the uncomment
        # the following if clause:
#        if wc.extension != "1":
#            print("You should use the extension '1' for the single experiment case")
#            sys.exit(-1)
    else:
        name = wc.sample+wc.extension
    return "sequences/{exp}/{sample}/{name}.{read}.allTrimmed.fq.gz".format(read=read, exp=wc.experiment, sample=wc.sample, name=name)
        

def getMergeInput(config, wc):
    return ["results/{exp}/{sample}/data/{s}.mapped.{mapper}.bam".format(exp=wc.experiment, sample=wc.sample, mapper=wc.mapper,s=s) for s in config["metafiles"][wc.sample]]


# will take a format string like "{experiment}/{sample}{extension}" and create absolute paths 
def format2path(fmt,df):
    headers = list(df)
    output = []
    for i in range(df.shape[0]):
        f = fmt
        for h in headers:
            sub = "\{"+h +"\}"
            rep = df[h][i]
            if rep != rep: rep='.X' #chec if nan
            f = re.sub(sub,rep,f)
        output.append(f)
    #make unique outputs
    output = list(set(output))
    return(output)



def read_sampleinfo2(conf):
    infile = conf["sampleinfo"]["file"]
    if not os.path.exists(infile):
        sys.exit("Error: File "+infile+" does not exist!")

    input = pd.read_csv(infile,dtype=str,na_values='')
    samples = input["sample"].unique()
    output = dict()
    # create one dict per sample
    for s in samples:
        output[s]=dict()
        output[s]["in"]=[]
        output[s]["out"]=[]
        subdf = input[input["sample"]==s]
        subdf.reset_index(inplace=True,drop=True)
        for i in range(df.shape[0]):        
            output[s]["experiment"]
            
    return(output)

def read_sampleinfo(conf):
    infile = conf["sampleinfo"]["file"]
    if not os.path.exists(infile):
        sys.exit("Error: File "+infile+" does not exist!")

    input = pd.read_csv(infile,dtype=str,na_values='')
    samples = input["sample"].unique()
    output = dict()
    # create one dict per sample
    for s in samples:
        output[s]=dict()
        subdf = input[input["sample"]==s]
        subdf.reset_index(inplace=True,drop=True)
        for fmt in ["outfmt","mergefmt"]:
            f = conf["sampleinfo"][fmt]
            output[s][fmt]=format2path(f,subdf)
    return(output)
