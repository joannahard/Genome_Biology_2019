#!/usr/bin/env python3

import pysam, sys

#class ReadsDb:
    

# returns a 
def createReadsDb(hetReadsFile, homReadsFile, Gfile, Sfile):
    retdb = []
    i = 0
    G = open(Gfile, "rt")
    for row in G:
        row = row.split()
        if row[0] == "CHROM":
            continue
        chr = row[0]
        Gsam = int(row[1])
        key = "{c}:{g}".format(c=chr, g=Gsam)
        retdb.append({"name": key, "chr": chr, "G": Gsam})
        i += 1
    G.close()

    S = open(Sfile, "rt")
    i = 0
    for row in S:
        row = row.split()
        if row[0] == "CHROM":
            continue
        chr = row[0]
        if chr != retdb[i]["chr"]:
            print("Mismatch on line {i}: {Gfile} chromosome {gchr} != {Sfile} chromosome {schr}\n".format(i=i, Gfile=Gfile, gchr = retdb[i]["chr"], Sfile=Sfile, schr=chr))
            sys.exit(-1)
        Ssam = int(row[1])
        key = "{prev}:{s}".format(prev=retdb[i]["name"], s=Ssam)
        retdb[i]["name"] = key
        retdb[i]["S"] = Ssam
    S.close()

    ali = {}
    ali["het"] = pysam.AlignmentFile(hetReadsFile, "rb")
    ali["hom"] = pysam.AlignmentFile(homReadsFile, "rb")

    # Notice that pysam works with 0-based pos, while inpout pos (originally from sam-format) are 1-based
    for post in retdb:
        for z in ["het","hom"]:
            post[z] = { "reads" :{}, "states" : {} }
            post[z]["reads"]["R"] = []
            post[z]["reads"]["A"] = []
            for read in iterCoveringReads(ali[z], chr = "{}".format(post["chr"]), G = post["G"]-1, S = post["S"]-1):
                allele, state = getState(read, post["G"]-1, post["S"]-1)
                if allele not in post[z]["states"].keys():
                    post[z]["states"][allele] = state
                elif state != post[z]["states"][allele]:
                    print("Error: mismatching states at {}: {} != {}, compared to previous state".format(allele, state, post["states"][allele]))
                    sys.exit(-1)
                post[z]["reads"][allele].append(read)
    return retdb
                

# Helper functions
#
# returns an iterator over reads covering positions G and S; Note! G and S are 0-based positions
def iterCoveringReads(ali, chr, G, S):
    for read in ali.fetch(contig = chr, start = G, end = S+1):
        if read.reference_start <= G and read.reference_end >= S:
            yield read

# returns a tuple with allele (R= ref ot A=alt) and the states of G and S; Note! G and S are 0-based positions
def getState(read, G, S):
    s = read.query_alignment_sequence
    t = dict([(t[1] , [t[0],t[2]]) for t in read.get_aligned_pairs(with_seq=True)])
    g = s[t[G][0]]
    s = s[t[S][0]]
    allele = "R" if g == t[G][1] else "A"
    return (allele, [g,s])

# returns a list with all polymorphic positions relative to the reference genome; ; Note! returned positions are 0-based positions
def getVariants(read):
    s = read.query_alignment_sequence
    t = dict([(t[1] , [t[0],t[2]]) for t in read.get_aligned_pairs(with_seq=True)])
    return [s[t[pos][0]] for pos in t.keys() if s[t[pos][0]]!=t[pos][1]]

# Dummy main
def main():
    db = createReadsDb("data/bam_files/fem1.bam", "data/bam_files/fib1.bam", "data/test_positions/G1.txt", "data/test_positions/S1.txt")
    print(db)

if __name__ == "__main__":
    main()

                                
# db = {}
# #for read in het.fetch(region="1:{G}-{S}".format(G=G_sam, S=S_sam)):
# for read in getCoveringReads(contig = "1", start = G, end = S):
#     key = getState(read, G)
#     if key not in db.keys():
#         db[key] = []
#     db[key].append(read)

# for g in db.keys():
#     for read in db[g]:
#         print("{G}{S}".format(S=getState(read,S),G=getState(read,G)))


