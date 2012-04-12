import bioLibCG
from cgNexus import Nexus
import cgDL
from cgAutoCast import autocast
from cgAutoKeyWord import autokey
from bioLibJA import subit
import cgWig

@autocast
def getReadsSpots(fN, chrom, strand, wigDir, outFN):

    #strand caste doesnt matter
    coord_eLevel = cgWig.loadSingleWigFloat(wigDir, chrom, strand, "ALL")
    
    gName_gInfo = {} #name: (numReads, numSpots)
    f = open(fN, 'r')
    fOut = open(outFN, 'w')
    for line in f:
        ls = line.strip().split('\t')
        geneName, gChrom, gStrand, gStarts, gEnds = ls

        gStarts = [int(x) for x in gStarts.split(',')]
        gEnds = [int(x) for x in gEnds.split(',')]
        exonPairs = zip(gStarts, gEnds)
        
        #get info
        totNumSpots = 0 #this is for L in 1/this
        totNumReads = 0 #NOTE: psuedo count is introduced in testing, not here
        totNumNSpots = 0 #this is for N for p val 
        for eStart, eEnd in exonPairs:
            for i in range(eStart, eEnd + 1):
                numReadsHere = coord_eLevel.get(i, 0)
                if numReadsHere >= 1:
                    totNumSpots += 1
                    totNumReads += numReadsHere

                if numReadsHere >= 3:
                    totNumNSpots += 1


        #write info
        fOut.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (geneName, gChrom, gStrand, totNumReads, totNumSpots, totNumNSpots) )


    fOut.close()
    f.close()

def getTotalSpots(allFN, formatFN):

    NX = Nexus(allFN, formatFN)
    NX.load(['numNSpots', 'numReads'])

    totalSpots = 0
    totalReads = 0
    while NX.nextID():

        totalSpots += NX.numNSpots
        totalReads += NX.numReads

    print 'spotsNSpots', totalSpots
    print 'numReads', totalReads



if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

