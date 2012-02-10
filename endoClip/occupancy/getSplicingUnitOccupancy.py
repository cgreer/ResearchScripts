import bioLibCG
import cgNexusFlat
import compareData
import cgWig
import matplotlib.pyplot as plt
import random
import numpy as np
import scipy.stats as stats
import scipy
from cgAutoCast import autocast

@autocast
def simulatePVal(totalUnitLength, numOnes, numTwos, actualOverlap, cutOffVal):
    '''cutOffVal is just for output for para'''

    #get sampling distribution for
    samplingDist = []
    for i in xrange(200): 
        ones = np.random.random_integers(1, totalUnitLength, numOnes)
        twos = np.random.random_integers(1, totalUnitLength, numTwos)
        ovlp = np.intersect1d(ones, twos)
        samplingDist.append(ovlp.size)

   
    sDist = np.array(samplingDist)
    theStd, theMean = sDist.std(), sDist.mean()
    zScore = (float(actualOverlap) - theMean)/theStd
    thePVal = scipy.special.ndtr(-zScore)

    print '\t'.join([str(x) for x in (theStd, theMean, zScore, thePVal, cutOffVal)])  
    
def getSplicingUnitOccupancy(tranFN, wigDir1, wigDir2, chrom, strand, maxCut):
        '''get the number of spots in each data set, and the number that overlap'''
        '''wigDir2 has to be hela cuz strand flip'''
        maxCut = int(maxCut)

        oppStrand = bioLibCG.switchStrand(strand)
        coord_value1 = cgWig.loadSingleWig(wigDir1, chrom, strand, 'ALL')
        coord_value2 = cgWig.loadSingleWig(wigDir2, chrom, oppStrand, 'ALL')

        # 0, 0, 0 = num1, num2, numOverlap
        covered = set()
        cutoff_overlap = dict( (i, [0, 0, 0]) for i in range(maxCut))
        f = open(tranFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                tChrom, tStrand = ls[1], bioLibCG.switchStrandFormat(ls[2])
                if tChrom != chrom or tStrand != strand:
                        continue
                tStart, tEnd = int(ls[3]), int(ls[4]) - 1
                cStart, cEnd = int(ls[5]), int(ls[6]) - 1
                exonStarts = [int(x) for x in ls[8][:-1].split(',')]
                exonEnds = [int(x) - 1 for x in ls[9][:-1].split(',')]
                exonPairs = zip(exonStarts, exonEnds)
                codingStatus = '_coding' in ls[13]
                tID = ls[0]

                #calulate intron pairs
                intronPairs = []
                i = 0
                for pair in exonPairs:
                        if i == 0:
                                i += 1
                                continue
                        iStart = exonPairs[i -1][1] + 1
                        iEnd = exonPairs[i][0] - 1
                        intronPairs.append((iStart, iEnd))
                        i += 1


                
                #take care of messy UTRs and assign utr ranges
                #5UTR
                if strand == '1':
                        if cStart == tStart or cStart == tEnd + 1:
                                range5 = ()
                        else:
                                range5 = (tStart, cStart - 1)
                else:
                        if cEnd + 1 == tStart or cEnd + 1 == tEnd + 1:
                                range5 = ()
                        else:
                                range5 = (cEnd + 1, tEnd)

                
                #3UTR
                if strand == '1':
                        if cEnd + 1 == tStart or cEnd + 1 == tEnd + 1:
                                range3 = ()
                        else:
                                range3 = (cEnd + 1, tEnd)
                else:
                        if cStart == tStart or cStart == tEnd + 1:
                                range3 = ()
                        else:
                                range3 = (tStart, cStart - 1)
                                
                utr5 = compareData.subtractTwoRanges([range5], intronPairs)
                utr3 = compareData.subtractTwoRanges([range3], intronPairs)
                
                exonPairs = compareData.subtractTwoRanges(exonPairs, [range5])
                exonPairs = compareData.subtractTwoRanges(exonPairs, [range3])

                pairs__type = [ (exonPairs, 'C_EXON'), (intronPairs, 'C_INTRON') ]
                for pairs, type in pairs__type:
                    for pair in pairs:
                        for i in xrange(pair[0], pair[1] + 1):
                                if codingStatus:
                                    if type == 'C_EXON':
                                        if i in covered: continue #multiple transcripts will have same exons
                                        covered.add(i)
                                        val1 = coord_value1.get(i, 0)
                                        val2 = coord_value2.get(i, 0)

                                        for cut in range(1, maxCut):
                                            #in1 = (val1 >= cut)
                                            #in2 = (val2 >= cut)
                                            in1 = (val1 == cut)
                                            in2 = (val2 == cut)

                                            if in1 and in2:
                                                cutoff_overlap[cut][2] += 1
                                            
                                            if in1:
                                                cutoff_overlap[cut][0] += 1

                                            if in2:
                                                cutoff_overlap[cut][1] += 1

                                    elif type == 'C_INTRON':
                                        #intronChr_strand_coord.setdefault(tChrom, {}).setdefault(tStrand, set()).add(i)
                                        pass

        for i in range(1, maxCut):

            cutoff_overlap[i].extend(['%s:%s' % (chrom, strand), i])
            pString = '\t'.join([ str(x) for x in cutoff_overlap[i] ])
            print pString

@autocast
def compactCounts(countsFN, outFN, cutMax):

    cut_oneTwoOverlap = dict( (i, [0,0,0]) for i in range(1,cutMax) )

    f = open(countsFN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        oneCount, twoCount, overlap, chromStrand, cut = ls
        oneCount, twoCount, overlap, cut = int(oneCount), int(twoCount), int(overlap), int(cut)

        cut_oneTwoOverlap[cut][0] += oneCount
        cut_oneTwoOverlap[cut][1] += twoCount
        cut_oneTwoOverlap[cut][2] += overlap
    f.close()

    pStrings = ['%s\t%s\t%s\t%s\n' % (i, cut_oneTwoOverlap[i][0], cut_oneTwoOverlap[i][1], cut_oneTwoOverlap[i][2]) for i in range(1, cutMax) ]
    f = open(outFN, 'w')
    f.writelines(pStrings)    
    f.close()
    
def plotCutOccupancy(plotInfo, imgName):

    ones = []
    twos = []
    overlaps = []
    f = open(plotInfo, 'r')
    for line in f:
        ls = line.strip().split('\t')
        cut, oneVal, twoVal, overlap = ls
        ones.append(int(oneVal))
        twos.append(int(twoVal))
        overlaps.append(int(overlap))
    f.close()
    

    plt.title('Occupancy of HeLa and U87')
    plt.ylabel('Number of Spots')
    plt.xlabel('Expression Cutoff for Spot')                             
    plt.plot(ones, color = 'r', label = 'U87')
    plt.plot(twos, color = 'b', label = 'HeLa')
    plt.plot(overlaps, color = 'k', label = 'Overlap')
    plt.legend()
    plt.show()



if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
