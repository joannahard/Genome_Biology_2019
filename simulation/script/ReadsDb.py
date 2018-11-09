#!/usr/bin/env python3

import pysam, sys, random, math

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
        for read in iterCoveringReads(hetReads, chr = self.chr, G = self.G-1, S = self.S-1):
            self.addRead("het", read)
        for read in iterCoveringReads(homReads, chr = self.chr, G = self.G-1, S = self.S-1):
            self.addRead("hom", read)

    def addRead(self, zyg, read):
        allele, state = getState(read, self.G - 1, self.S - 1)
        if state != self.zyg[zyg].states[allele] and self.zyg[zyg].states[allele] != None:
            print("Error: mismatching states at {}: {} != {}, compared to previous state".format(allele, state, self.zyg[zyg].states[allele]))
            sys.exit(-1)
        self.zyg[zyg].states[allele] = state
        self.zyg[zyg].reads[allele].append(read)

    def simulateLocus(self, T, SNV, EAL, fADO, locusCounts):
        C = unnest(T)
        zygReads = { "l" : { c : self.zyg["hom"] for c in C } }

        # sSNV or not
        if SNV:
            v = random.choice(T)
            for c in v:
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

        


class ReadsDb:
    def __init__(self, hetReadsFile, homReadsFile, Gfile, Sfile):
        self.loci = []
        self.fill(hetReadsFile, homReadsFile, Gfile, Sfile)

    def iter_simulate(self, T, pSNV, pEAL, f_ADO, locusCounts):
        L = len(self.loci)
        SNV = [(random.random() < pSNV) for l in range(L) ]  # S is SNV
        EAL = [(random.random() < pEAL) for l in range(L) ] # existence of alignment error
        for l in range(len(self.loci)):
            yield self.loci[l].simulateLocus(T, SNV[l], EAL[l], 0.5, locusCounts) 

    def simulate(self, T, pSNV, pEAL, f_ADO, locusCounts):
        L = len(self.loci)
        SNV = [(random.random() < pSNV) for l in range(L) ]  # S is SNV
        EAL = [(random.random() < pEAL) for l in range(L) ] # existence of alignment error
        return [ self.loci[l].simulateLocus(T, SNV[l], EAL[l], 0.5, locusCounts) for l in range(len(self.loci)) ]

    def fill(self, hetReadsFile, homReadsFile, Gfile, Sfile):
        i = 0
        G = open(Gfile, "rt")
        S = open(Sfile, "rt")
        for rowG, rowS in zip(G, S):
            rowG = rowG.split()
            if rowG[0] == "CHROM":
                continue
            rowS = rowS.split()
            chr = rowG[0]
            if chr != rowS[0]:
                print("non-matching chromosomes for line {i}: file {gf} -- {gc} != {sc} -- file {sf}".format(i=i, gc=chr, sc=rowS[0],
                                                                                                             gf= Gfile, sf=Sfile))
                sys.exit(-1)
            Gsam = int(rowG[1])
            Ssam = int(rowS[1])
            key = "{c}:{g}".format(c=chr, g=Gsam)
            self.loci.append(Locus(chr, Gsam, Ssam))
            i += 1
        G.close()
        S.close()
        hetReads = pysam.AlignmentFile(hetReadsFile, "rb")
        homReads = pysam.AlignmentFile(homReadsFile, "rb")

        # Notice that pysam works with 0-based pos, while inpout pos (originally from sam-format) are 1-based
        for locus in self.loci:
            locus.addAllReads(hetReads, homReads)
            
    



def main():
    ## Example usage
    T = { 0:{0,1}, 1:{2,3} }
    print(T)
    print(unnest(T))

    pSNV = 0.5
    pEAL = 0.5
    fADO = 0.5
    
    locusCounts = { "Ref": {1 : 10000, 2 : 0}, "Mut": {1 : 10000, 2 : 0} }
    readsDB = ReadsDb("data/bam_files/fem1.bam", "data/bam_files/fib1.bam", "data/test_positions/G1.txt", "data/test_positions/S1.txt")
    
    reads = readsDB.simulate(T, pSNV, pEAL, fADO, locusCounts)

    for c in unnest(T):
        print("cell {}".format(c))
        n=0
        for l in reads:
            print("locus {}, {} reads".format(n, len(l[c])))
            n+=1
            for r in l[c]:
                print(r)

        header = { 'HD': {'VN': '1.0'},
                   'SQ': [{'LN': 1575, 'SN': 'chr1'},
                      {'LN': 1584, 'SN': 'chr2'}] }
        
        myReads = pysam.AlignmentFile("myReads_{}.bam".format(c), "wb", header=header)
        for l in reads:
            for r in l[c]:
                myReads.write(r)
                           


            
# Helper functions of possible use
#
# returns an iterator over reads covering positions G and S; Note! G and S are 0-based positions
# returns a list with all polymorphic positions relative to the reference genome; ; Note! returned positions are 0-based positions
def getVariants(read):
    s = read.query_alignment_sequence
    t = dict([(t[1] , [t[0],t[2]]) for t in read.get_aligned_pairs(with_seq=True)])
    return [s[t[pos][0]] for pos in t.keys() if s[t[pos][0]]!=t[pos][1]]

def iterCoveringReads(reads, chr, G, S):
    for read in reads.fetch(contig = chr,
                            start = G,
                            end = S+1):
        if read.reference_start <= G and read.reference_end >= S:
            yield read

# returns a tuple with allele (R= ref ot A=alt) and the states of G and S; Note! G and S are 0-based positions
def getState(read, G, S):
    s = read.query_alignment_sequence
    t = dict([(t[1] , [t[0],t[2]]) for t in read.get_aligned_pairs(with_seq=True)])
    g = s[t[G][0]]
    s = s[t[S][0]]
    allele = "Ref" if g == t[G][1] else "Mut"
    return (allele, [g,s])
 
if __name__ == "__main__":
    main()

                                
