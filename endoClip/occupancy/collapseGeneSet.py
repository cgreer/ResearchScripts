import bioLibCG
import cgNexusFlat
import compareData
from interval import IntervalSet, Interval
from cgAutoCast import autocast

@autocast
def createCollapsedGeneSets(tranFN, outFN, acceptableTypes = 'EXON', onlyCoding = True):
    '''get areas occupied by all transcripts in a gene'''

    acceptableTypes = acceptableTypes.strip().split(',')

    geneName_intervalSet = {}
    geneName_info = {}
    f = open(tranFN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        tChrom, tStrand = ls[1], bioLibCG.switchStrandFormat(ls[2])
        #if tChrom != chrom or tStrand != strand:
                #continue
        tStart, tEnd = int(ls[3]), int(ls[4]) - 1
        cStart, cEnd = int(ls[5]), int(ls[6]) - 1
        exonStarts = [int(x) for x in ls[8][:-1].split(',')]
        exonEnds = [int(x) - 1 for x in ls[9][:-1].split(',')]
        exonPairs = zip(exonStarts, exonEnds)
        codingStatus = '_coding' in ls[13]
        geneName = ls[10]
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
        if tStrand == '1':
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
        if tStrand == '1':
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

        geneName_info.setdefault(geneName, set()).add((tChrom, tStrand))
        pairs__type = [ (exonPairs, 'EXON'), (intronPairs, 'INTRON'), (utr5, '5UTR'), (utr3, '3UTR') ]
        for pairs, type in pairs__type:
            for pair in pairs:
                if type in acceptableTypes:
                    
                    if onlyCoding and not codingStatus: continue

                    #create geneset/info if does not exist
                    if geneName not in geneName_intervalSet:
                        geneName_intervalSet[geneName] = IntervalSet()


                    geneName_intervalSet[geneName].add(Interval(pair[0], pair[1] + 1))

    for geneName, info in geneName_info.iteritems():
        if len(info) > 1:
            if geneName in geneName_intervalSet:
                del geneName_intervalSet[geneName] # if it spans different chromosomes/strands...

    fOut = open(outFN, 'w')
    for geneName, iSet in geneName_intervalSet.iteritems():
        gStarts = []
        gEnds = []
        
        for interv in iSet:
            gStarts.append(interv.lower_bound)
            gEnds.append(interv.upper_bound)

        chrom, strand = geneName_info[geneName].pop()
        outString = [geneName, chrom, strand, ','.join([str(x) for x in gStarts]), ','.join([str(x) for x in gEnds])]             
        fOut.write('\t'.join([str(x) for x in outString]) + '\n')
            


    fOut.close()

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
