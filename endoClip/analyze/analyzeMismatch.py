import bioLibCG
import cgNexusFlat
import cgOriginRNAFlat
import cgAlignmentFlat
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats



def mismatchPositionBoxPlot(oFN, aFN, minSNR = 2):

        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['filteredTargets', 'snrSS'])
        
        aNX = cgNexusFlat.Nexus(aFN, cgAlignmentFlat.cgAlignment)
        aNX.load(['mismatchPositions'])

        position_mCount = {}
        for oID in oNX.ids:

                #filter
                if oNX.snrSS[oID] < float(minSNR):                        
                        continue
                
                for aID in oNX.filteredTargets[oID]:

                        for mPosition in aNX.mismatchPositions[aID]:
                                position_mCount[mPosition] = position_mCount.get(mPosition, 0) + 1

        counts = []
        for pos in sorted(position_mCount.keys()):
                print pos, position_mCount[pos]
                counts.append([position_mCount[pos]])

        plt.boxplot(counts)
        plt.show()

        
if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
