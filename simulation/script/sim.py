#!/usr/bin/env python3


import pysam, sys, random, math
from Locus import *



def unnest(adict):
    if type(adict) in [ int, str, float ]:
        return [adict]
    ret = []
    if type(adict) == dict:
        for i in adict.values():
            ret.extend(unnest(i))
    else:
        for i in adict:
            ret.extend(unnest(i))        
    return ret
        
    

def simulateLocus(locusReads, T, SNV, EAL, fADO, locusCounts):
    if type(locusReads) != Locus:
        print("TypeError: argument locusReads must be of type Locus")
        sys.exit(-1)

    C = unnest(T)
    zygReads = { "l" : { c : locusReads.zyg["hom"] for c in C } }

    # sSNV or not
    if SNV:
        v = random.choice(T)
        for c in v:
            zygReads["l"][c] = locusReads.zyg["het"]
            
    fReads = { "l" : { c : { "Ref" :  1 , "Mut" : 1 } for c in C } }

    # Aligment error
    if EAL:
        if SNV:
            zygReads["lE"] = { c : locusReads.zyg["hom"] for c in C }
        else:
            zygReads["lE"] = { c : locusReads.zyg["het"] for c in C }

        fReads["lE"] = { c : { "Ref" : 1, "Mut" : 1 } for c in C }


    #ADO on l, randomly on either allele
    for c in random.sample(C, math.floor(len(C) * fADO)):
        allele = "Ref" if random.random() > 0.5 else "Mut"
        fReads["l"][c][allele] = 0

    # #ADO on lE, randomly on either allele
    # for i in random.sample(locusCounts.keys(), fADO * locusCOunts.keys()):
    #     allele = "Ref" if random.random > 0.5 else "Mut"
    #     freads["lE"][c][allele] = 0

    # Allelic imbalance, modeled by locusCounts
    for zyg in fReads.keys():
        for c in C:
            for allele in ["Ref", "Mut"]:
                counts = list(locusCounts[allele].keys())
                weights = list(locusCounts[allele].values())
                fReads[zyg][c][allele] *= random.choices(counts, weights=weights, k=1)[0]

    # output
    ret = { c : [] for c in C }
    for c in C:
        for zyg in fReads.keys():
            ret[c].extend(zygReads[zyg][c].sample(fReads[zyg][c]["Ref"],fReads[zyg][c]["Mut"]))

    return ret



def main():
    T = { 0:{0,1}, 1:{2,3} }
    print(T)
    print(unnest(T))

    pSNV = 0.5
    pEAL = 0.5
    
    readsDB = createReadsDb("data/bam_files/fem1.bam", "data/bam_files/fib1.bam", "data/test_positions/G1.txt", "data/test_positions/S1.txt")
    L= len(readsDB)
    SNV = [(random.random() < pSNV) for l in range(L) ]  # S is SNV
    EAL = [(random.random() < pEAL) for l in range(L) ] # existence of alignment error
    fADO = 0.5
    locusCounts = { "Ref": {1 : 10000, 2 : 0}, "Mut": {1 : 10000, 2 : 0} }
    reads = [ simulateLocus(readsDB[l], T, SNV[l], EAL[l], 0.5, locusCounts) for l in range(L) ]

    for c in unnest(T):
        print("cell {}".format(c))
        n=0
        for l in reads:
            print("locus {}, {} reads".format(n, len(l[c])))
            n+=1
            # for r in l[c]:
            #     print(r)
        
 
if __name__ == "__main__":
    main()

