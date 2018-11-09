#!/usr/bin/env python3

import pysam, sys, random


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

    def addRead(self, zyg, read):
        allele, state = getState(read, self.G - 1, self.S - 1)
        if state != self.zyg[zyg].states[allele] and self.zyg[zyg].states[allele] != None:
            print("Error: mismatching states at {}: {} != {}, compared to previous state".format(allele, state, self.zyg[zyg].states[allele]))
            sys.exit(-1)
        self.zyg[zyg].states[allele] = state
        self.zyg[zyg].reads[allele].append(read)

    def addAllReads(self,  hetReads, homReads):
        for read in iterCoveringReads(hetReads, chr = self.chr, G = self.G-1, S = self.S-1):
            self.addRead("het", read)
        for read in iterCoveringReads(homReads, chr = self.chr, G = self.G-1, S = self.S-1):
            self.addRead("hom", read)

    def __str__(self):
        ret = "" + \
        "chr = {}\n".format(self.chr) + \
        "G = {}\n".format(self.G) + \
        "S = {}\n".format(self.S) + \
        ""
        return ret
        



# returns a list of Locus instances
def createReadsDb(hetReadsFile, homReadsFile, Gfile, Sfile):
    ret = []
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
            print("non-matching chromosomes for line {i}: file {gf} -- {gc} != {sc} -- file {sf}".format(i=i, gc=chr, sc=rowS[0], gf= Gfile, sf=Sfile))
            sys.exit(-1)
        Gsam = int(rowG[1])
        Ssam = int(rowS[1])
        key = "{c}:{g}".format(c=chr, g=Gsam)
        ret.append(Locus(chr, Gsam, Ssam))
        i += 1
    G.close()
    S.close()
    hetReads = pysam.AlignmentFile(hetReadsFile, "rb")
    homReads = pysam.AlignmentFile(homReadsFile, "rb")
    
    # Notice that pysam works with 0-based pos, while inpout pos (originally from sam-format) are 1-based
    for locus in ret:
        locus.addAllReads(hetReads, homReads)
    return ret
                

# Helper functions
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


# Dummy main
def main():
    db = createReadsDb("data/bam_files/fem1.bam", "data/bam_files/fib1.bam", "data/test_positions/G1.txt", "data/test_positions/S1.txt")
    for locus in db:
        print(locus)

if __name__ == "__main__":
    main()

                                
