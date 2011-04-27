import bioLibCG
import cgPlot
import cgPeaks

def makeFigure(fN, targetFN, alignmentFN, cName):
        #make targetDict
        f = open(targetFN, 'r')
        targetDict = {} # tID: tLoc
        for line in f:
                ls = line.strip().split('\t')
                targetDict[int(ls[0])] = ls[1]
        f.close()

        #make alignmentDict
        alignDict = {} # sid: {target: offset}
        f = open(alignmentFN, 'r')
        for line in f:
                ls = line.strip().split(' ')
                sID = int(ls[0])
                tID = int(ls[1])
                offset = int(ls[4])
                if not sID in alignDict:
                           alignDict[sID] = {}

                alignDict[sID][tID] = offset #assumes one source to target...
        f.close()

        f = open(fN, 'r')
        
        histoVals = []

        for line in f:
                ls = line.strip().split('\t')
                sID = int(ls[0])
                sLoc = ls[1]
                sChrom, sStrand, sStart, sEnd = bioLibCG.tccSplit(sLoc)
                sLen = sEnd - sStart
                tIDs = ls[4].split(',')

                for tID in tIDs:
                        tID = int(tID)
                        tLoc = targetDict[tID]
                        chrom, strand, start, end = bioLibCG.tccSplit(tLoc)
                        offset = alignDict[sID][tID]

                        if sStrand == '1':
                                start = start - 19 + offset
                                end = start + sLen
                        else:
                                end = end + 19 - offset
                                start = end - sLen

                        scanRange = bioLibCG.makeTcc(chrom, strand, start, end)
                                
                        stretch = cgPeaks.stretch(scanRange, cName)
                        highest = stretch.getHighestLevel()
                        sortedKeys = stretch.profile.keys()

                        if sStrand == '-1':
                                sortedKeys.reverse()

                        i = 0
                        for key in sortedKeys:
                                level = stretch.profile[key]
                                for j in range(0,level):
                                        histoVals.append(i)
                                i += 1




        cgPlot.plotHistogram(histoVals)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(makeFigure, sys.argv)
