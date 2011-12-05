import bioLibCG
import cgNexusFlat
import cgOriginRNAFlat
import cgAlignmentFlat
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats
import math

def eLevelHistogram(oFN, aFN, oRNA = True):
        
        oRNA = 'True' in oRNA

        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['filteredTargets', 'eLevel'])

        aNX = cgNexusFlat.Nexus(aFN, cgAlignmentFlat.cgAlignment)
        aNX.load(['tELevel'])
     
        histValues = []
        for oID in oNX.eLevel:
                if oRNA:
                        histValues.append(oNX.eLevel[oID])
                else:
                        for aID in oNX.filteredTargets[oID]:
                                histValues.append(aNX.tELevel[aID])

        histVals = [math.log(x, 10) for x in histValues]
        plt.hist(histVals, 50)
        type = 'oRNA'
        if not oRNA: type = 'Targets (degradome)'
        plt.title('Expression Level for %s' % type)
        plt.xlabel('log(Expression Level)')
        plt.ylabel('Number of %s' % type)


        plt.show()



if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
