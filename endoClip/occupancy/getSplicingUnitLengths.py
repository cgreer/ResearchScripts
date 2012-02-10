import bioLibCG
import cgNexusFlat
import compareData
from cgAutoCast import autocast

def getSplicingUnitLengths(tranFN, wigDir, chrom, strand):

        
        exonChr_strand_coord = {}
        intronChr_strand_coord = {}
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
                                        exonChr_strand_coord.setdefault(tChrom, {}).setdefault(tStrand, set()).add(i)
                                    elif type == 'C_INTRON':
                                        #intronChr_strand_coord.setdefault(tChrom, {}).setdefault(tStrand, set()).add(i)
                                        pass

        iLength = 0
        for chrom in intronChr_strand_coord:
            for strand in intronChr_strand_coord[chrom]:
                iLength += len(intronChr_strand_coord[chrom][strand])

        eLength = 0
        for chrom in exonChr_strand_coord:
            for strand in exonChr_strand_coord[chrom]:
                eLength += len(exonChr_strand_coord[chrom][strand])

        print 'total Exon Length (all exons overlapped)', eLength                        
        print 'total intron Length (all introns overlapped)', iLength                        





if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
