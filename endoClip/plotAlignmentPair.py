import bioLibCG
import cgAlignmentFlat
import cgOriginRNAFlat
import matplotlib.pyplot as plt
import cgWig
import cgNexusFlat
from matplotlib.ticker import MaxNLocator

def plotPair(oValues, aValues, maxO, maxA, tStart, tEnd):
        
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(oValues, 'b-')

        #plot the small NTs
        smallNT = [oValues[x] if tStart <= x <= tEnd else -10 for x in range(len(oValues))]
        ax1.plot(smallNT, 'bo', markersize = 9)
        ax1.set_ylim(0, 70)
        ax1.yaxis.set_major_locator(MaxNLocator(7))
        ax1.set_yticklabels( [str((float(x)/20) * maxO) for x in [0, 10, 20, 30, 40, 50, 60, 70]] )
        ax1.set_xlabel('nt of small/target')
        ax1.set_ylabel('Expression Level (origin RNA)')
        ax1.set_title('oRNA/Target Expression Pair')

        ax2 = ax1.twinx()
        ax2.plot(aValues, 'r-')
        ax2.set_ylim(-30, 70)
        ax2.set_yticklabels( [str((float(x)/50) * maxA) for x in [-30,-20, -10, 0, 10, 20, 30, 40, 50, 60, 70]] )
        ax2.set_ylabel('Expression Level (target RNA)')
        ax2.yaxis.set_major_locator(MaxNLocator(10))
        
        ax2.xaxis.set_major_locator(MaxNLocator(len(oValues) - 1))
        ax1.xaxis.set_major_locator(MaxNLocator(len(oValues) - 1))
        print len(oValues)
        ax2.set_xlim(0, len(oValues))

        
        plt.show()

def plotByIDs(oID, aID, plotData):

        f = open(plotData, 'r')
        for line in f:
                ls = line.strip().split('\t')
                oIDF = ls[0]
                aIDF = ls[1]
                oValues = [float(x) for x in ls[2].split(',')] 
                aValues = [float(x) for x in ls[3].split(',')] 
                maxO = int(ls[4]) 
                maxA = int(ls[5])
                tStart = int(ls[6])
                tEnd = int(ls[7])

                if oID == oIDF and aID == aIDF:
                        plotPair(oValues, aValues, maxO, maxA, tStart, tEnd)
        f.close() 


def getPairInfo(oFN, aFN, oWigDir, aWigDir, outFN):

        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['tcc', 'filteredTargets'])

        aNX = cgNexusFlat.Nexus(aFN, cgAlignmentFlat.cgAlignment)
        aNX.load(['tTcc', 'tStart', 'tEnd', 'tLength' ])
        
        #get expression Dicts for alignments and small RNAs               
        oWigDict = cgWig.loadWigDict(oWigDir)
        aWigDict = cgWig.loadWigDict(aWigDir)
        
        fOut = open(outFN, 'w')
        for oID in oNX.tcc:
                for aID in oNX.filteredTargets[oID]:


                        #expand the oTcc to fit the target
                        oChrom, oStrand, oStart, oEnd = bioLibCG.tccSplit(oNX.tcc[oID])
                        oStart -= aNX.tStart[aID]
                        oEnd += (aNX.tLength[aID] - 1) - aNX.tEnd[aID]
                        oTcc = bioLibCG.makeTcc(oChrom, oStrand, oStart, oEnd)

                        #expand the peak tcc to be the target tcc (peak is only 1-4nt long)
                        aChrom, aStrand, aStart, aEnd = bioLibCG.tccSplit(aNX.tTcc[aID])
                        aStart -= 25
                        aEnd += 25
                        aTcc = bioLibCG.makeTcc(aChrom, aStrand, aStart, aEnd)
                        
                        print oNX.tcc[oID], oTcc, aNX.tTcc[aID], aTcc
                        
                        #get expression
                        oCoord_value = cgWig.getExpressionProfile(oTcc, oWigDict)
                        aCoord_value = cgWig.getExpressionProfile(aTcc, aWigDict)
                        
                        oValues = [oCoord_value[x] for x in sorted(oCoord_value.keys())]
                        aValues = [aCoord_value[x] for x in sorted(aCoord_value.keys())]

                        maxO = max(oValues[aNX.tStart[aID]:aNX.tEnd[aID]]) # might want to change to max only within the small range
                        maxA = max(aValues)
                      
                        print oValues
                        print aValues
                        oValues = [(float(x)/maxO) * 20 for x in oValues]
                        aValues = [(float(x)/maxA) * 50 for x in aValues]

                        pString = [str(oID), str(aID), ','.join([str(x) for x in oValues]), ','.join([str(x) for x in aValues]), str(maxO), str(maxA), str(aNX.tStart[aID]), str(aNX.tEnd[aID])]
                        fOut.write('\t'.join(pString) + '\n')

        fOut.close()                        


        

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
