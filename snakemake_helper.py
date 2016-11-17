#!/usr/bin/env python3

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
