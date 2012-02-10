import bioLibCG
import cgNexusFlat
import cgWig
from cgAutoCast import autocast
import scipy.stats.stats as pStats
import numpy as np
import random

def fullSpanFromPairs(listOfPairs):

    span = []
    for pairStart, pairEnd in listOfPairs:
        span.extend(range(pairStart, pairEnd))

    return span        

def mixSpanByBin(theSpan, binSize):
    '''bins will be non-overlapping'''

    maxSpanIndex = len(theSpan)
    binCoords = list(np.arange(0, maxSpanIndex, binSize))
    if maxSpanIndex not in binCoords: binCoords.append(maxSpanIndex)
    start__end = zip(binCoords[:-1], binCoords[1:])
    random.shuffle(start__end)

    newSpan = []
    for start, end in start__end:
        newSpan.extend(theSpan[start:end])

    return newSpan        

def binAvg(coord_val, coords):

    theSum = sum([coord_val.get(x, 0) for x in coords])
    avg = float(theSum)/len(coords)
    return avg

@autocast
def byGene(geneSpanFN, wigDir1, wigDir2, chrom, strand, outFN, simulation = False):
    '''hela must be 2nd wigDir2 cuz strand flip'''
    strand = str(strand) #undo autocast
   
    print 'loading wigs'
    oppStrand = bioLibCG.switchStrand(strand)
    coord_value1 = cgWig.loadSingleWig(wigDir1, chrom, strand, 'ALL')
    coord_value2 = cgWig.loadSingleWig(wigDir2, chrom, oppStrand, 'ALL')

    print 'calculating bin values'
    f = open(geneSpanFN, 'r')
    fOut = open(outFN, 'w')
    for line in f:
        ls = line.strip().split('\t')
        sChrom, sStrand = ls[1], ls[2]
        if sChrom != chrom or sStrand != strand:
            continue
        geneName = ls[0]
        geneStarts = [int(x) for x in ls[3].split(',')]
        geneEnds = [int(x) for x in ls[4].split(',')]
        spanPairs = zip(geneStarts, geneEnds)

        frameLength = 10
        skipAmount = 2
        theSpan = fullSpanFromPairs(spanPairs)
        spanLength = len(theSpan)


        binAvgs1 = []
        binAvgs2 = []

        for theBinAvg, theCoord_Val in [(binAvgs1, coord_value1), (binAvgs2, coord_value2)]:
            #mix up bins if simulation
            if simulation:
                newSpan = mixSpanByBin(theSpan, frameLength)
            else:
                newSpan = theSpan
            
            i = 0
            while (i+frameLength) < (spanLength+1):
                binNums = newSpan[i:(i + frameLength)]
                theBinAvg.append(binAvg(theCoord_Val, binNums))
                i = i + skipAmount

        #get rid of all 0,0 pairs for correlation 
        editPairs = zip(binAvgs1, binAvgs2)
        newPairs = [pair for pair in editPairs if not (pair[0] == 0 and pair[1] == 0)]
        newX = [pair[0] for pair in newPairs]
        newY = [pair[1] for pair in newPairs]

        dataLoad = sum(binAvgs1) + sum(binAvgs2)
        dataLoad = float(dataLoad)/2
        pcc = pStats.pearsonr(binAvgs1, binAvgs2)
        scc, pVal = pStats.spearmanr(binAvgs1, binAvgs2)
        outString = [geneName, pcc[0], ','.join([str(x) for x in binAvgs1]), ','.join([str(x) for x in binAvgs2]), '%s:%s:%s' % (sChrom, sStrand, theSpan[0]), dataLoad, scc]  
        fOut.write('\t'.join([str(x) for x in outString]) + '\n')

    fOut.close()
    f.close()
    

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

