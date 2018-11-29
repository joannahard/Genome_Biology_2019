#!/usr/bin/env python3

## This file implements the simlation model/algorithm described in ../doc/sim.text/pdf


import pysam, sys, random, math, json
from collections import Counter, OrderedDict


class Reads:
    def __init__(self):
        self.reads = { "Ref" : [],
                  "Mut" : [] }
        self.states = { "Ref" : None,
                   "Mut" : None }
    def sample(self, nRef, nMut):
        try:
            return random.choices(self.reads["Ref"], k=nRef) + random.choices(self.reads["Mut"], k=nMut)
        except IndexError:
            sys.stdout.write("Error: No reads to sample from\n")
            sys.stdout.write("{m} Mut reads and "
                             "{r} Ref reads available\n".format(m=len(self.reads["Mut"]),
                                                                r=len(self.reads["Ref"])))
            sys.exit(-1)

    def __iter__(self):#iterAllReads():
        try:
            for read in self.reads["Ref"] + self.reads["Mut"]:
                yield read
        except IndexError:
            sys.stdout.write("Error: No reads to iterate over\n")
            sys.exit(-1)

    def asList(self):#iterAllReads():
        return  self.reads["Ref"] + self.reads["Mut"]

        
class Locus:
    def __init__(self, chr, G, S, flanksize):
        self.chr = "{}".format(chr)
        self.G = int(G)
        self.S = int(S)
        self.zyg = { "het" : Reads(),
                     "hom" : Reads() }
        self.flanksize = flanksize
        self.flanks = { "het" : [],
                        "hom" : [] }


    def fullName(self):
        ret = "{}:{}-{}".format(self.chr, self.G, self.S)
        # ret = "" + \
        # "chr = {}\n".format(self.chr) + \
        # "G = {}\n".format(self.G) + \
        # "S = {}\n".format(self.S) + \
        # ""
        return ret

    def __str__(self):
        ret = "{}:{}".format(self.chr, self.S)
        # ret = "" + \
        # "chr = {}\n".format(self.chr) + \
        # "G = {}\n".format(self.G) + \
        # "S = {}\n".format(self.S) + \
        # ""
        return ret

    def setStates(self, zyg, readset):
        obs_states = {"Ref": [], "Mut" : []}
        for read in iterCoveringReads(readset, chr = self.chr,
                                      G = self.G-1, S = self.S-1):
            if read.reference_start <= self.G-1 and read.reference_end >= self.S:
                allele, state = getState(read, self.G-1, self.S-1)
                obs_states[allele].append(state)
        # We need to check for occasional aberrant states (seq errors?)
        for allele in ["Ref", "Mut"]:
            tmp = Counter(obs_states[allele])
            if len(tmp) > 1:
                sys.stdout.write("Abberrant '{a}' states found at "
                                 "{chr}:{G}-{S}, for the '{z}' zygote; "
                                 "the most common is chosen as the true "
                                 "'{a}' state:\n{o}\n".format(a = allele,
                                                            chr=self.chr,
                                                            G=self.G,
                                                            S=self.S,
                                                            z=zyg, o=tmp))
            self.zyg[zyg].states[allele] = max(tmp,key=tmp.get)
        
    def addAllReads(self, hetReads, homReads):
        self.setStates("het", hetReads)
        n = 0
        for read in iterCoveringReads(hetReads, chr = self.chr,
                                      G = self.G-self.flanksize-1,
                                      S = self.S+self.flanksize-1):
            self.addRead("het", read,n)
            n+=1
        self.setStates("hom", homReads)
        if n==0:
            sys.exit(-1)
        n = 0
        for read in iterCoveringReads(homReads, chr = self.chr,
                                      G = self.G-self.flanksize-1,
                                      S = self.S+self.flanksize-1):
            self.addRead("hom", read,n)
            n+=1
        if n==0:
            sys.exit(-1)

    def addRead(self, zyg, read, n):
        if read.reference_start <= self.G-1 and read.reference_end >= self.S:
            allele, state = getState(read, self.G - 1, self.S - 1)
            if state in self.zyg[zyg].states[allele]:
                self.zyg[zyg].reads[allele].append(read)
            else:
                sys.stdout.write("Error: Abberrant state found for "
                                 "read {n} at {chr}:{G}-{S}, '{z}' "
                                 "allele '{a}': {s} != {p}, compared "
                                 "to previous state. This read is "
                                 "ignored\n".format(n=n, chr=self.chr,
                                                  G=self.G, S=self.S,
                                                  z=zyg, a=allele, s=state,
                                                  p=self.zyg[zyg].states[allele]))
        else:
            self.flanks[zyg].append(read)

    def simulateLocus(self, T, SNV, EAL, fADO, locusCounts, ADOfreqs = True):
        C = unnest(T)
        zygReads = { "l" : { c : self.zyg["hom"] for c in C } }
        flanks = { c : self.flanks["hom"] for c in C }

        # sSNV or not
        if SNV:
            L= flatten(T)
            v = random.choice(list(L.keys()))
            for c in L[v]:
                zygReads["l"][c] = self.zyg["het"]
                flanks[c] = self.flanks["het"]
                

        fReads = { "l" : { c : { "Ref" :  1 , "Mut" : 1 } for c in C } }

        # Aligment error
        if EAL:
            if SNV:
                zygReads["lE"] = { c : self.zyg["hom"] for c in C }
            else:
                zygReads["lE"] = { c : self.zyg["het"] for c in C }

            fReads["lE"] = { c : { "Ref" : 1, "Mut" : 1 } for c in C }


        #ADO on l, randomly on either allele
        # if ADOfreqs = true then sample that a fraction fADO of cells, else
        # use fADO as iid probability of each cell having ADO
        adoC = [ [] for c in C ]
        if ADOfreqs:
            sample = random.sample(C, math.floor(len(C) * fADO))
            for c in sample:
                zyg = "l" if random.random() < 0.5 else "lE"
                allele = "Ref" if random.random() < 0.5 else "Mut"
                fReads["l"][c][allele] = 0
                adoC[c].append("{}:{}".format(zyg,allele))
        else:
            for c in C:  
                for zyg in fReads.keys():
                    for allele in [ "Ref", "Mut"]:
                        if random.random() < fADO:
                            fReads[zyg][c][allele] = 0
                            adoC[c].append("{}:{}".format(zyg,allele))


        # #ADO on lE, randomly on either allele
        # for i in random.sample(locusCounts.keys(), fADO * locusCounts.keys()):
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
        states = {}
        reads = { c : [] for c in C }
        for c in C:
            states["cell{}".format(c)] = zygReads["l"][c].states
            for zyg in fReads.keys():
                reads[c].extend(zygReads[zyg][c].sample(fReads[zyg][c]["Ref"],fReads[zyg][c]["Mut"]))
                

        return { 'locus' : self.fullName(), 'SNV' : SNV, 'EAL' : EAL,
                 'ADO' : adoC,
                 #{ "cell{}".format(c): True if c in adoC else False for c in C },
                 'states' : states, 'reads' : reads, 'flanks' : flanks }

    def allReadsAndFlanks(self):
        return self.zyg["hom"].asList() + self.flanks["hom"]

        


class ReadsDb:
    def __init__(self, hetReadsFile, homReadsFile, Sitesfile, flanksize = 1000):
        self.loci = []
        self.header = None
        self.fill(hetReadsFile, homReadsFile, Sitesfile, flanksize)

    def simulate(self, T, f_SNV, f_EAL, f_ADO, locusCounts, freqs = True):
        L = range(len(self.loci))
        SNV = EAL = None
        # if freqs = true then sample that fraction fSNV/fEAL of loci, else
        # use fSNV/fEAL as iid probabilities of each locus being SNV/EAL
        if freqs:
            # existence of alignment error
            EAL, SNV = assignSubSample(L, f_EAL, f_SNV) 
        else:
            # existence of alignment error
            EAL = [ (random.random() < f_EAL) for l in L ] 
            SNV = [ False if EAL[l] else (random.random() < f_SNV) for l in L ]  # S is SNV only if not having EAL
            
        return { str(self.loci[l]) :
                 self.loci[l].simulateLocus(T, SNV[l], EAL[l],
                                            f_ADO, locusCounts,
                                            freqs)
                 for l in range(len(self.loci)) }

    
    def simulateAndWriteToFile(self, T, f_SNV, f_EAL, f_ADO, locusCounts,
                               outprefix, flanks = True, freqs=False):
        
        reads = self.simulate(T, f_SNV, f_EAL, f_ADO, locusCounts, freqs)
        with open("{o}_genVals.json".format(o=outprefix), "wt") as genVals:
            vals = { k : { i:v[i] for i in ['locus', 'SNV', 'EAL',
                                           'ADO', 'states'] }
                                                     for k,v in reads.items() }
            genVals.write("{}\n".format(json.dumps(vals, indent=2)))
        
        for c in unnest(T):
            rgtag = "cell{}".format(c)
            myheader = self.header
            for tag in ['ID','LB','SM']:
                myheader['RG'][0][tag]= rgtag
                
            alfile = "{o}_cell{c}.bam".format(o=outprefix, c=c)
            cellReads = pysam.AlignmentFile(alfile, "wb", header=myheader)
            for l in reads.values():
                for r in l['reads'][c]:
                    r.set_tag("RG",rgtag)
                    cellReads.write(r)
                if flanks == True:
                    for r in l['flanks'][c]:
                        r.set_tag("RG",rgtag)
                        cellReads.write(r)


    def writeBulkToFile(self, outprefix):
        myheader = self.header
        for tag in ['ID','LB','SM']:
            myheader['RG'][0][tag]= "bulk"
        myReads = pysam.AlignmentFile("{o}_bulk.bam".format(o=outprefix), "wb", header=myheader)
        for locus in self.loci:
            for read in locus.allReadsAndFlanks():
                read.set_tag("RG","bulk")
                myReads.write(read)
        
    def fill(self, hetReadsFile, homReadsFile, Sitesfile, flanksize):
        with open(Sitesfile, "rt") as sites:
            for row in sites:
                row = row.strip().split()
                if row[0] == "chrom":
                    continue
                chr = row[0]
                Gsam = int(row[1])
                Ssam = int(row[2])
                key = "{c}:{g}".format(c=chr, g=Gsam)
                self.loci.append(Locus(chr, Gsam, Ssam, flanksize))
        hetReads = pysam.AlignmentFile(hetReadsFile, "rb")
        self.header=hetReads.header.as_dict()
        homReads = pysam.AlignmentFile(homReadsFile, "rb")
        
        # Notice that pysam works with 0-based pos, while input pos
        # (originally from sam-format) are 1-based
        for locus in self.loci:
            locus.addAllReads(hetReads, homReads)
            
    



def main():
    # TODO: Make use of argsparse modelue instead
    if len(sys.argv) != 9:
        print("Usage: {} <Reads-file Heterozygotes> <Reads-file Homozygotes> <sites-file> <coverage-file> <counts_limit> <popTree> <f_SNV> <f_EAL> <f_ADO> <output reads file>)".format(sys.argv[0]))
        sys.exit(-1)

    hetfile = sys.argv[1]
    homfile = sys.argv[2]
    sitesfile = sys.argv[3]
    count_file = sys.argv[4]
    counts_limit = sys.argv[5]
    popT = sys.argv[6]
    f_SNV = sys.argv[7]
    f_EAL = sys.argv[8]
    f_ADO = sys.argv[9]
    outprefix = sys.argv[10]

    locusCounts = createLocusCounts(counts_file, counts_limit)
    
    readsDB = ReadsDb(hetfile, homfile, sitesfile)
    
    readsDB.simulateAndWriteToFile(popT, f_SNV, f_EAL, f_ADO, LocusCounts, outprefix)

    readsDB.writeBulkToFile(outprefix)

    
            
# Helper functions 
#
def unnest(adict):
    if type(adict) in [ int, str, float ]:
        return [adict]
    ret = []
    if type(adict) in [dict, OrderedDict]:
        for i in adict.values():
            ret.extend(unnest(i))
    else:
        for i in adict:
            ret.extend(unnest(i))        
    return ret

def flatten(adict):
    ret = {}
    if type(adict) in [dict, OrderedDict]:
        for i in adict.keys():
            ret.update(flatten(adict[i]))
            ret[i] = unnest(adict[i])
    return ret

def assignSubSample(L, f_EAL, f_SNV):
    k = math.floor(f_EAL *len(L))
    sample = random.sample(L, k)
    EAL = [ True if l in sample else False for l in L ]
    
    tmp = [ l for l in L if l not in sample ]
    k = math.floor(f_SNV *len(tmp))
    sample = random.sample(tmp, k)
    SNV = [ True if l in sample else False for l in L ]
    
    return (EAL, SNV)

# def assignSubSample(L, f):
#     k = math.floor(f *len(L))
#     sample = random.sample(L, k)
#     return [ True if l in sample else False for l in L ]

def createLocusCounts(counts_file, limit):
    locusCounts = { "Ref": [], "Mut": [] }
    head = {}
    nsites=0
    with open(counts_file, 'rt') as counts:
        for row in counts:
            row = row.strip().split()
            if len(head) == 0:
                for i,val in enumerate(row):
                    head[val] = i
                nsites = int((len(head)-3)/5)
            else:
                dp = sum([ int(row[head["{}_DP".format(i)]]) > 0  for i in range(1,nsites+1) ])
                if dp  > 12:
                    for a,b in [ ("Ref","REF"), ("Mut","ALT") ]:
                        nt = row[head[b]]
                        locusCounts[a].extend([ int(row[head["{i}_{nt}".format(i=1,nt=nt)]]) for i in range(1,nsites+1) if int(row[head["{i}_{nt}".format(i=1,nt=nt)]]) > limit ])
    locusCounts = { key : Counter(locusCounts[key]) for key in locusCounts.keys() }
    return locusCounts

def iterCoveringReads(reads, chr, G, S):
    if G < 0:
        print("G={}".format(G))
        G = 0
    if S >= reads.get_reference_length(chr):
        print("S={}".format(S))
        S = reads.get_reference_length(chr)-1
    for read in reads.fetch(contig = chr,
                            start = G,
                            end = S+1):
        yield read

# returns an iterator over reads covering positions G and S; Note! G and S are 0-based positions
# returns a list with all polymorphic positions relative to the reference genome; ; Note! returned positions are 0-based positions
def getVariants(read):
    s = read.query_alignment_sequence
    t = dict([(t[1] , [t[0],t[2]]) for t in read.get_aligned_pairs(with_seq=True)])
    return { pos:(s[t[pos][0]],t[pos][1]) for pos in t.keys() if s[t[pos][0]]!=t[pos][1] }


# returns a tuple with allele (R= ref ot A=alt) and the states of G and S; Note! G and S are 0-based positions
def getState(read, G, S):
    s = read.query_alignment_sequence
    sstart = read.query_alignment_start
    t = dict([(t[1], [t[0],t[2]]) for t in read.get_aligned_pairs(with_seq=True)])
    g = s[t[G][0] - sstart]
    s = s[t[S][0] - sstart]
    allele = "Ref" if g == t[G][1] else "Mut"
    return (allele, "_".join([g,s]))


if __name__ == "__main__":
    main()

                                
