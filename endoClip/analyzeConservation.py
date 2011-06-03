import bioLibCG
import cgNexusFlat
import cgOriginRNAFlat
import cgAlignmentFlat
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats

def conservedHisto(oFN):
        
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['phastScores'])

        scores = []
        for oID in oNX.phastScores:

                avgScore = sum(oNX.phastScores[oID]) / float(len(oNX.phastScores[oID]))

                               
                #scores.extend(oNX.phastScores[oID])
                scores.append(avgScore)

        print len(scores)
        plt.title('PhastCons Scores By Nucleotide')
        plt.title('PhastCons Scores By oRNA')
        
        plt.ylabel('Number of Nucleotides')
        plt.ylabel('Number of oRNA')
        
        plt.xlabel('PhastCons Score')
        plt.xlabel('PhastCons Score Average')
        
        plt.hist(scores, 50)
        plt.show()

def correlationSNR(oFN):
        
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['phastScores', 'snrSS'])

        snrX = []
        scoreY = []
        for oID in oNX.phastScores:
                
                snr = oNX.snrSS[oID]
                avgScore = sum(oNX.phastScores[oID]) / float(len(oNX.phastScores[oID]))
                
                snrX.append(snr)
                scoreY.append(avgScore)
                '''                
                for pScore in oNX.phastScores[oID]:
                        snrX.append(snr)
                        scoreY.append(pScore)
                '''            

        conserved = [True if x > .8 else False for x in scoreY]
        print conserved.count(True)
        print len(conserved)
        plt.title('SNR vs PhastCons Score')
        
        plt.ylabel('Avg PhastCons Score of oRNA')
        
        plt.xlabel('SNR')
        
        plt.plot(snrX, scoreY, 'ro')
        plt.show()

def phastScoreByNT(oFN, aFN, oIDFilter = None):

        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['phastScores', 'snrSS', 'filteredTargets'])
        
        aNX = cgNexusFlat.Nexus(aFN, cgAlignmentFlat.cgAlignment)
        aNX.load(['mismatchPositions'])

        misPositionsX = []
        misScoresY = []
        positionsX = []
        scoresY = []

        for oID in oNX.phastScores:

                avgScore = sum(oNX.phastScores[oID]) / float(len(oNX.phastScores[oID]))
                
                #filter
                if (avgScore < .90) or (oNX.snrSS[oID] < 2):
                        continue
               
                if oIDFilter:
                        if oID != int(oIDFilter):
                                continue

                misPositions = set()
                #get consolidated mismatches
                for aID in oNX.filteredTargets[oID]:
                        for mPos in aNX.mismatchPositions[aID]:

                                misPositions.add(mPos)

                for i, pScore in enumerate(oNX.phastScores[oID]):
                        if i in misPositions:
                                misPositionsX.append(i + 1.2) #1 is for 0BASE, .2 is for differentiating between mis and reg
                                misScoresY.append(pScore)
                        else:
                                positionsX.append(i + 1)
                                scoresY.append(pScore)
        
        highestNT = max(positionsX)
        for i in range(1, highestNT + 1):
                plt.axvspan(i - .15, i + .35, facecolor='g', alpha=.25)

        plt.plot(positionsX, scoresY, 'bo')
        plt.plot(misPositionsX, misScoresY, 'ro')
        plt.ylim(0, 1.1)
        plt.xlim(0, 24)
        plt.title('Conservation by Position of Conserved oRNA')
        plt.ylabel('PhastCons Score')
        plt.xlabel('Nucleotide Position')

        
        plt.show()

def boxplotConservation(oFN, minAvgPhast = .90, minSNR = 2, oRNAType = 'oRNA'):

        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['phastScores', 'snrSS'])
        
        position_pScores = {}
        for oID in oNX.phastScores:

                avgScore = sum(oNX.phastScores[oID]) / float(len(oNX.phastScores[oID]))
                
                #filter
                if (avgScore < float(minAvgPhast)) or (oNX.snrSS[oID] < float(minSNR)):
                        continue
                
                if avgScore == 1.00: continue
                for i, pScore in enumerate(oNX.phastScores[oID]):
                             position_pScores.setdefault(i, []).append(pScore)           

        data = []
        for pos in sorted(position_pScores.keys()):
                data.append(position_pScores[pos])

        plt.boxplot(data)
        plt.ylim(.9, 1.1)
        highestPosition = max(position_pScores.keys()) 
        plt.xlim(0, highestPosition + 2)
        plt.title('Conservation by Position of Conserved %s (avg cons > %s)' % (oRNAType, minAvgPhast) )
        plt.ylabel('PhastCons Score')
        plt.xlabel('Nucleotide Position')
       
        for i in range(len(data)):
                plt.text(i + 1, 1.05, str(len(data[i])))
        
        plt.show()

def testNTDifference(oFN, minAvgPhast = .90, minSNR = 2, oRNAType = 'oRNA'):

        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['phastScores', 'snrSS'])
        
        groupA = [10,11,12,13]
        #groupB = [15,16,17,18]
        groupB = [4,5,6,7]
        
        a = []
        b = []

        for oID in oNX.phastScores:

                avgScore = sum(oNX.phastScores[oID]) / float(len(oNX.phastScores[oID]))
                
                #filter
                if (avgScore < float(minAvgPhast)) or (oNX.snrSS[oID] < float(minSNR)):
                        continue
                
                if avgScore == 1.00: continue
                for i, pScore in enumerate(oNX.phastScores[oID]):
                        if (i + 1) in groupA:
                                a.append(pScore)

                        if (i + 1) in groupB:
                                b.append(pScore)
        
        print len(a)/4, len(b)/4
        print stats.ranksums(a,b)                                



if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
