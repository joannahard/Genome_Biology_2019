#!/usr/bin/env python3

import pysam, sys, random, math
from collections import Counter

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
            print("Error: No reads to sample from")
            sys.exit(-1)
            

class Locus:
    def __init__(self, chr, G, S):
        self.chr = "{}".format(chr)
        self.G = int(G)
        self.S = int(S)
        self.zyg = { "het" : Reads(),
                "hom" : Reads() }

    def __str__(self):
        ret = "" + \
        "chr = {}\n".format(self.chr) + \
        "G = {}\n".format(self.G) + \
        "S = {}\n".format(self.S) + \
        ""
        return ret

    def addAllReads(self,  hetReads, homReads):
        n = 0
        for read in iterCoveringReads(hetReads, chr = self.chr,
                                      G = self.G-1, S = self.S-1):
            self.addRead("het", read,n)
            n+=1
        n = 0
        for read in iterCoveringReads(homReads, chr = self.chr,
                                      G = self.G-1, S = self.S-1):
            self.addRead("hom", read,n)
            n+=1

    def addRead(self, zyg, read, n):
        allele, state = getState(read, self.G - 1, self.S - 1)
        if state != self.zyg[zyg].states[allele] and self.zyg[zyg].states[allele] != None:
            print("Error: mismatching states for read {n} at {chr}:{G}-{S}, {z} allele {a}: {s} != {p}, compared to previous state".format(n=n, chr=self.chr, G=self.G, S=self.S, z=zyg, a=allele, s=state, p=self.zyg[zyg].states[allele]))
            return
            sys.exit(-1)
        self.zyg[zyg].states[allele] = state
        self.zyg[zyg].reads[allele].append(read)

    def simulateLocus(self, T, SNV, EAL, fADO, locusCounts):
        C = unnest(T)
        zygReads = { "l" : { c : self.zyg["hom"] for c in C } }

        # sSNV or not
        if SNV:
            v = random.choice(list(T.keys()))
            for c in T[v]:
                zygReads["l"][c] = self.zyg["het"]

        fReads = { "l" : { c : { "Ref" :  1 , "Mut" : 1 } for c in C } }

        # Aligment error
        if EAL:
            if SNV:
                zygReads["lE"] = { c : self.zyg["hom"] for c in C }
            else:
                zygReads["lE"] = { c : self.zyg["het"] for c in C }

            fReads["lE"] = { c : { "Ref" : 1, "Mut" : 1 } for c in C }


        #ADO on l, randomly on either allele
        for c in random.sample(C, math.floor(len(C) * fADO)):
            allele = "Ref" if random.random() > 0.5 else "Mut"
            fReads["l"][c][allele] = 0

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
        ret = { c : [] for c in C }
        for c in C:
            for zyg in fReads.keys():
                ret[c].extend(zygReads[zyg][c].sample(fReads[zyg][c]["Ref"],fReads[zyg][c]["Mut"]))

        return ret

        


class ReadsDb:
    def __init__(self, hetReadsFile, homReadsFile, Sitesfile):
        self.loci = []
        self.template = None
        self.fill(hetReadsFile, homReadsFile, Sitesfile)

    def iter_simulate(self, T, pSNV, pEAL, f_ADO, locusCounts):
        L = range(len(self.loci))
        SNV = assignSubSample(L, f_SNV) # S is SNV
        EAL = assignSubSample(L, f_EAL) # existence of alignment error
        # L = len(self.loci)
        # SNV = [(random.random() < pSNV) for l in range(L) ]  # S is SNV
        # EAL = [(random.random() < pEAL) for l in range(L) ] # existence of alignment error
        for l in range(len(self.loci)):
            yield self.loci[l].simulateLocus(T, SNV[l], EAL[l], 0.5, locusCounts) 

    def simulate(self, T, f_SNV, f_EAL, f_ADO, locusCounts):
        L = range(len(self.loci))
        SNV = assignSubSample(L, f_SNV) # S is SNV
        EAL = assignSubSample(L, f_EAL) # existence of alignment error
        # L = len(self.loci)
        # SNV = [(random.random() < pSNV) for l in range(L) ]  # S is SNV
        # EAL = [(random.random() < pEAL) for l in range(L) ] # existence of alignment error
        return [ self.loci[l].simulateLocus(T, SNV[l], EAL[l], 0.5, locusCounts) for l in range(len(self.loci)) ]

    def simulateAndWriteToFile(self, T, f_SNV, f_EAL, f_ADO, locusCounts, outprefix):
        reads = self.simulate(T, f_SNV, f_EAL, f_ADO, locusCounts)
        header = { 'HD': {'VN': '1.0'},
                   'SQ': [{'LN': 1575, 'SN': 'chr1'},
                          {'LN': 1584, 'SN': 'chr2'}] }

        for c in unnest(T):
            myReads = pysam.AlignmentFile("{o}_cell{c}.bam".format(o=outprefix, c=c), "wb", template=self.template)
            for l in reads:
                for r in l[c]:
                    myReads.write(r)

    
    def fill(self, hetReadsFile, homReadsFile, Sitesfile):
        with open(Sitesfile, "rt") as sites:
            for row in sites:
                row = row.strip().split()
                if row[0] == "chrom":
                    continue
                chr = row[0]
                Gsam = int(row[1])
                Ssam = int(row[2])
                key = "{c}:{g}".format(c=chr, g=Gsam)
                self.loci.append(Locus(chr, Gsam, Ssam))
        hetReads = pysam.AlignmentFile(hetReadsFile, "rb")
        self.template=hetReads
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


            
# Helper functions 
#
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

def assignSubSample(L, f):
    k = math.floor(f *len(L))
    sample = random.sample(L, k)
    return [ 1 if l in sample else 0 for l in L ]

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
    for read in reads.fetch(contig = chr,
                            start = G,
                            end = S+1):
        if read.reference_start <= G and read.reference_end >= S+1:
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
    return (allele, [g,s])
 
if __name__ == "__main__":
    main()

                                
