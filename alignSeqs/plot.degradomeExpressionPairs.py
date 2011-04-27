import bioLibCG
import cgDB
import cgAlignment
import cgOriginRNA
import cgPeaks
import matplotlib.pyplot as plt

def plotPairs(oDir, aDir, cName):

        oDC = cgDB.dataController(oDir, cgOriginRNA.OriginRNA)
        id_oRNA = oDC.load()

        aDC = cgDB.dataController(aDir, cgAlignment.cgAlignment)
        id_alignment = aDC.load()

        for oID, oRNA in id_oRNA.items():
                
                if not oRNA.passedFilter:
                        continue

                for aID in oRNA.filteredTargets:

                        alignment = id_alignment[aID]
                        chrom, strand, start, end = bioLibCG.tccSplit(alignment.tTcc)
                        offset = alignment.tStart
                        sLen = alignment.sLength
                        print sLen
                        print oRNA.sequence
                        print oRNA.tcc
                        print alignment.tTcc
                        if strand == '1':
                                start = start - 19 + offset
                                end = start + sLen
                        else:
                                end = end + 19 - offset
                                start = end - sLen

                        print chrom, strand, start, end
                        scanRange = bioLibCG.makeTcc(chrom, strand, start, end)
                        
                        stretch = cgPeaks.stretch(scanRange, cName)
                        sortedKeys = stretch.profile.keys()
                        sortedKeys.sort()

                        if strand == '-1':
                                sortedKeys.reverse()
                        

                        xVals = range(1, sLen + 2)
                        xVals = sortedKeys
                        yVals = [stretch.profile[x] for x in sortedKeys]
                        print xVals, len(xVals)
                        print yVals, len(yVals)
                        
                        plt.plot(xVals, yVals)
                        plt.show()

                        return 0

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(plotPairs, sys.argv)

