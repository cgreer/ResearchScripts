import bioLibCG
import cgAlignmentFlat
import cgOriginRNAFlat
import matplotlib.pyplot as plt
import cgWig
import cgNexusFlat
from matplotlib.ticker import MaxNLocator
import GenomeFetch as gf
import numpy as np

def plotPair2(oValues, aValues, tStart, tEnd, oTcc, aTcc, oSeq, dSeq, misPositions, oContext, aContext, imgName):
        #everything should be set up like the alignment...oRNA is 5 to 3'

        maxO = max(oValues[tStart:tEnd + 1]) # might want to change to max only within the small range
        maxA = max(aValues)
      
        oValues = [(float(x)/maxO) * 20 for x in oValues]
        aValues = [(float(x)/maxA) * 50 for x in aValues]
     
     
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
      
        #plot the small NTs
        smallNT = [oValues[x] if tStart <= x <= tEnd else -10 for x in range(len(oValues))]
        ax1.plot(smallNT, 'ko', markersize = 9)
        
        #plot the small expression line + letters
        ax1.plot(oValues, 'b-')
        
        #t is x, s is y ---> fill 
        t = range(len(oValues))
        s = oValues
        ax1.fill_between(t, -10, s, color = 'b')
        
        #plot small letters
        for i, val in enumerate(oValues):
                plt.text(i, val, oSeq[i])
      
        #small conservation
        #phastScores = [1.0 for i in range(tStart, tEnd + 1)]
        #phastScores[0] = 0.0
        #for i, a in enumerate(range(tStart, tEnd + 1)):
                #plt.text(a, 10 + 2*phastScores[i] - .4, '*')
                #plt.text(a, 10, '--')
                #plt.text(a, 12, '--')

        #plot sequence above
        plt.text(tStart - 1, 66, '5\'')
        plt.text(tEnd + 1, 66, '3\'')
        plt.text(tStart - 1, 68, '3\'')
        plt.text(tEnd + 1, 68, '5\'')
        print misPositions
        for i in range(tStart, tEnd + 1):
                plt.text(i, 66, oSeq[i], color = 'b')
                plt.text(i, 68, bioLibCG.compSeq(dSeq[i]), color = 'r')
                if i not in misPositions:
                        plt.axvline(x = i + .12, ymin = .96, ymax = .97, color='k')
                else:
                        plt.text(i, 67, 'x')

        ax1.set_ylim(0, 70)
        ax1.yaxis.set_major_locator(MaxNLocator(7))
        ax1.set_yticklabels( [str((float(x)/20) * maxO) for x in [0, 10, 20, 30, 40, 50, 60, 70]] )
        ax1.set_xlabel('nt of small/target')
        ax1.set_ylabel('Expression Level (origin RNA)')
        ax1.set_title('oRNA/Target Expression Pair')

        #plot deg line, letters
        ax2 = ax1.twinx()
        dLine = ax2.plot(aValues, 'k-')
        
        
        #t is x, s is y ---> fill 
        t = range(len(aValues))
        s = aValues
        ax2.fill_between(t, 0, s, color = 'r')

        for i, val in enumerate(aValues):
                plt.text(i, val, dSeq[i])
        ax2.set_ylim(-30, 70)
        ax2.set_yticklabels( [str((float(x)/50) * maxA) for x in [-30,-20, -10, 0, 10, 20, 30, 40, 50, 60, 70]] )
        ax2.set_ylabel('Expression Level (target RNA)')
        ax2.yaxis.set_major_locator(MaxNLocator(10))

       
        #x axis stuff
        ax2.xaxis.set_major_locator(MaxNLocator(len(oValues) - 1))
        ax1.xaxis.set_major_locator(MaxNLocator(len(oValues) - 1))
        ax2.set_xlim(0, len(oValues))
        
        #expression boundaries
        #cStart = tStart + 9 - 2
        #cEnd = tStart + 10 + 2
        #plt.axvline(x = cStart, color='g')
        #plt.axvline(x = cEnd, color = 'g')
        #plt.axvspan(cStart, cEnd, facecolor='g', alpha=.25)

        #info box
        #plt.text(42, 60, 'tcc: %s' % oTcc, color = 'b')
        #plt.text(42, 58, 'gene: %s' % 'None', color = 'b')
        #plt.text(42, 56, 'context: %s' % oContext, color = 'b')
        #plt.text(42, 54, 'type: %s' % 'None', color = 'b')

        #plt.text(42, 50, 'tcc: %s' % aTcc, color = 'r')
        #plt.text(42, 48, 'gene: %s' % 'None', color = 'r')
        #plt.text(42, 46, 'context: %s' % aContext, color = 'r')
        #plt.text(42, 44, 'type: %s' % 'None', color = 'r')
        
        #plt.savefig(imgName, dpi=300, bbox_inches='tight', pad_inches=.22)
        plt.show()
        plt.close('all')


def plotByIDs(oID, aID, plotData, imageDir):

        f = open(plotData, 'r')
        for line in f:
                ls = line.strip().split('\t')
                oIDF = ls[0]
                aIDF = ls[1]
                oValues = [float(x) for x in ls[2].split(',')] 
                aValues = [float(x) for x in ls[3].split(',')] 
                tStart = int(ls[4])
                tEnd = int(ls[5])
                oTcc = ls[6]
                aTcc = ls[7]
                oSeq = ls[8]
                dSeq = ls[9]
                oSNR = ls[11]
                misPositions = []
                oContext = ls[12]
                aContext = ls[13]
                if ls[10] != '':
                        misPositions = [tStart + int(x) for x in ls[10].split(',')]

                if oID == oIDF and aID == aIDF:
                        imageName = imageDir + '/%s.%s.%s.png' % (oSNR, oID, aID)
                        plotPair2(oValues, aValues, tStart, tEnd, oTcc, aTcc, oSeq, dSeq, misPositions, oContext, aContext, imageName)
        f.close() 

def plotAllIDs(sortedFN, imageDir):

    f = open(sortedFN, 'r')
    for line in f:
            ls = line.strip().split('\t')
            oIDF = ls[0]
            aIDF = ls[1]
            oValues = [float(x) for x in ls[2].split(',')] 
            aValues = [float(x) for x in ls[3].split(',')] 
            tStart = int(ls[4])
            tEnd = int(ls[5])
            oTcc = ls[6]
            aTcc = ls[7]
            oSeq = ls[8]
            dSeq = ls[9]
            oSNR = ls[11]
            misPositions = []
            oContext = ls[12]
            aContext = ls[13]
            if ls[10] != '':
                    misPositions = [tStart + int(x) for x in ls[10].split(',')]

            imageName = imageDir + '/%s.%s.%s.png' % (oSNR, oIDF, aIDF)
            plotPair2(oValues, aValues, tStart, tEnd, oTcc, aTcc, oSeq, dSeq, misPositions, oContext, aContext, imageName)
    f.close() 

def getPairInfo2(oFN, aFN, oWigDir, aWigDir, outFN, assembly):

        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['snrSS', 'tcc', 'filteredTargets', 'context'])

        aNX = cgNexusFlat.Nexus(aFN, cgAlignmentFlat.cgAlignment)
        aNX.load(['tTcc', 'tStart', 'tEnd', 'tLength', 'mismatchPositions', 'context'])
        
        #get expression Dicts for alignments and small RNAs               
        oWigDict = cgWig.loadWigDict(oWigDir)
        aWigDict = cgWig.loadWigDict(aWigDir)
       
        myG = gf.GenomeFetch(assembly)

        fOut = open(outFN, 'w')
        for oID in oNX.tcc:
                for aID in oNX.filteredTargets[oID]:


                        #expand the oTcc to fit the target
                        oChrom, oStrand, oStart, oEnd = bioLibCG.tccSplit(oNX.tcc[oID])
                        if oStrand == '1':
                                oStart -= aNX.tStart[aID] #5' 
                                oEnd += (aNX.tLength[aID] - 1) - aNX.tEnd[aID] #3'
                                oTcc = bioLibCG.makeTcc(oChrom, oStrand, oStart, oEnd)
                        elif oStrand == '-1':                                
                                oStart -= (aNX.tLength[aID] - 1) - aNX.tEnd[aID] #3'
                                oEnd += aNX.tStart[aID] #5'
                                oTcc = bioLibCG.makeTcc(oChrom, oStrand, oStart, oEnd)
                                
                        #expand the peak tcc to be the target tcc (peak is only 1-4nt long)
                        aChrom, aStrand, aStart, aEnd = bioLibCG.tccSplit(aNX.tTcc[aID])
                        aStart -= 25
                        aEnd += 25
                        aTcc = bioLibCG.makeTcc(aChrom, aStrand, aStart, aEnd)
                        
                        
                        #get expression
                        oCoord_value = cgWig.getExpressionProfile(oTcc, oWigDict)
                        aCoord_value = cgWig.getExpressionProfile(aTcc, aWigDict)
                        
                        oValues = [oCoord_value[x] for x in sorted(oCoord_value.keys())]
                        aValues = [aCoord_value[x] for x in sorted(aCoord_value.keys())]

                        if oStrand == '-1':
                                oValues.reverse()
                        if aStrand == '-1':
                                aValues.reverse()

                        #get sequences
                        aSeq = myG.getSequence(aTcc)
                        oSeq = myG.getSequence(oTcc)

                        oSNR = oNX.snrSS[oID]
                        oContext = oNX.context[oID]
                        aContext = aNX.context[aID]
                        
                        pString = [str(oID),
                                   str(aID),
                                   ','.join([str(x) for x in oValues]),
                                   ','.join([str(x) for x in aValues]),
                                   str(aNX.tStart[aID]),
                                   str(aNX.tEnd[aID]),
                                   oNX.tcc[oID],
                                   aNX.tTcc[aID],
                                   oSeq,
                                   aSeq,
                                   ','.join([str(x) for x in aNX.mismatchPositions[aID]]),
                                   str(oSNR),
                                   oContext,
                                   aContext]
                       
                        fOut.write('\t'.join(pString) + '\n')

        fOut.close()                        

        

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
