from interval import IntervalSet, Interval
import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast
from geneProperty import geneProperty
import GenomeFetch as gf


def getGeneLength(geneRanges, genePropFN):
   
    NX = cgNexusFlat.Nexus(genePropFN, geneProperty)
    NX.load(['geneName', 'geneLength'])
   
    #make inverse dictionary
    gName_nID = {}
    for id in NX.ids:
        gName_nID[NX.geneName[id]] = id

    f = open(geneRanges, 'r')
    for line in f:
        ls = line.strip().split('\t')
        
        sChrom, sStrand = ls[1], ls[2]
        geneName = ls[0]
        geneStarts = [int(x) for x in ls[3].split(',')]
        geneEnds = [int(x) for x in ls[4].split(',')]
        spanPairs = zip(geneStarts, geneEnds)
        totalLength = sum([pair[1] - pair[0] for pair in spanPairs])

        #get id
        nID = gName_nID.get(geneName, None)

        if nID:
            NX.geneLength[nID] = totalLength
        
    f.close()
    
    NX.save()
    
def getGeneGC(genePropFN, rn = None, tn = None):
   
    myGF = gf.GenomeFetch('hg19')

    NX = cgNexusFlat.Nexus(genePropFN, geneProperty)
    NX.load(['geneName', 'geneChrom', 'geneStrand', 'geneStarts', 'geneEnds', 'geneGCContent'], [rn, tn])
   
    for id in NX.ids:

        spanPairs = zip(NX.geneStarts[id], NX.geneEnds[id])
        spanTccs = ['%s:%s:%s:%s' % (NX.geneChrom[id], NX.geneStrand[id], pair[0], pair[1]) for pair in spanPairs]

        totalSequenceLength = 0
        numGC = 0
        for tcc in spanTccs:
            seq = myGF.getSequence(tcc)
            totalSequenceLength += len(seq)
            numGC += seq.count('G') + seq.count('C')

        GCC = float(numGC)/totalSequenceLength

        NX.geneGCContent[id] = GCC
        

    NX.save()

def updateGeneData(geneRanges, genePropFN):
   
    myGF = gf.GenomeFetch('hg19')

    NX = cgNexusFlat.Nexus(genePropFN, geneProperty)
    NX.load(['geneName', 'geneChrom', 'geneStrand', 'geneStarts', 'geneEnds'])
   
    #make inverse dictionary
    gName_nID = {}
    for id in NX.ids:
        gName_nID[NX.geneName[id]] = id

    f = open(geneRanges, 'r')
    for line in f:
        ls = line.strip().split('\t')
        
        sChrom, sStrand = ls[1], ls[2]
        geneName = ls[0]
        geneStarts = [int(x) for x in ls[3].split(',')]
        geneEnds = [int(x) for x in ls[4].split(',')]

        #get id
        nID = gName_nID.get(geneName, None)

        if nID:
            NX.geneChrom[nID] = sChrom
            NX.geneStrand[nID] = sStrand
            NX.geneStarts[nID] = geneStarts
            NX.geneEnds[nID] = geneEnds
        
    f.close()

    NX.save()

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])


