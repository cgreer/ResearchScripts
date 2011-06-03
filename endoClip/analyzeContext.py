import bioLibCG
import cgNexusFlat
import cgOriginRNAFlat
import cgAlignmentFlat
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats

def pieFractions(countList):
        
        s = sum(countList)
        fracs = [float(x)/s for x in countList]
        return fracs

def positionsMappedHist(oFN, rn = None, tn = None):
        '''for results with high snrs, how many other places did this result map?'''

        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['tccs'], [rn, tn])

        histVals = []
        for oID in oNX.tccs:
                histVals.append(len(oNX.tccs[oID]))

        plt.title('Number of Mapping Locations')
        plt.ylabel('Number of oRNA')
        plt.xlabel('Number of Places Mapped')                             
        plt.hist(histVals, 20)
        plt.show()

def oRNAContextPie(oFN, rn = None, tn = None):
        '''REMEMBER!!! Have to do with grouped results...''' 
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['context2'], [rn, tn])

        context_count = {}

        for oID in oNX.context2:

                con = oNX.context2[oID]
                print con
                context_count[con] = context_count.get(con, 0) + 1

        labels = sorted(context_count.keys())
        counts = [context_count[x] for x in labels]
        fracs = pieFractions(counts)

        #add numbers to labels
        labels = ['%s (%s)' % (x, context_count[x]) for x in labels]
        
        #plot
        plt.pie(fracs, labels=labels, shadow = True)
        plt.show()

def oRNATypePie(oFN, rn = None, tn = None):        
        
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['transcriptType'], [rn, tn])

        context_count = {}

        for oID in oNX.transcriptType:

                con = oNX.transcriptType[oID]
                print con
                context_count[con] = context_count.get(con, 0) + 1

        labels = sorted(context_count.keys())
        counts = [context_count[x] for x in labels]
        fracs = pieFractions(counts)

        #add numbers to labels
        labels = ['%s (%s)' % (x, context_count[x]) for x in labels]
        
        #plot
        plt.pie(fracs, labels=labels, shadow = True)
        plt.show()

def targetContextPie(oFN, aFN, oContext = None, oType = None, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['context2'], [rn, tn])

        aNX = cgNexusFlat.Nexus(aFN, cgAlignmentFlat.cgAlignment)
        aNX.load(['context'], [rn, tn])

        context_count = {}

        for oID in oNX.context2:

                #filter oID types here
                con = oNX.context2[oID]
                typ = oNX.transcriptType[oID]
                
                if oContext:
                        if oContext != con: continue

                if oType:
                        if oType != typ: continue

                
                #gather targets' context info if oRNA is okay
                for aID in oNX.filteredTargets[oID]:
                        aCon = aNX.context[aID]
                        context_count[aCon] = context_count.get(aCon, 0) + 1

        labels = sorted(context_count.keys())
        counts = [context_count[x] for x in labels]
        fracs = pieFractions(counts)

        #add numbers to labels
        labels = ['%s (%s)' % (x, context_count[x]) for x in labels]
        
        #plot
        plt.pie(fracs, labels=labels, shadow = True)
        plt.show()

def targetTypePie(oFN, aFN, oContext = None, oType = None, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['transcriptType', 'context'], [rn, tn])

        aNX = cgNexusFlat.Nexus(aFN, cgAlignmentFlat.cgAlignment)
        aNX.load(['type'], [rn, tn])

        context_count = {}

        for oID in oNX.context:

                #filter oID types here
                con = oNX.context[oID]
                typ = oNX.transcriptType[oID]
                
                if oContext:
                        if oContext != con: continue

                if oType:
                        if oType != typ: continue

                
                #gather targets' context info if oRNA is okay
                for aID in oNX.filteredTargets[oID]:
                        aCon = aNX.type[aID]
                        context_count[aCon] = context_count.get(aCon, 0) + 1

        labels = sorted(context_count.keys())
        counts = [context_count[x] for x in labels]
        fracs = pieFractions(counts)

        #add numbers to labels
        labels = ['%s (%s)' % (x, context_count[x]) for x in labels]
        
        #plot
        plt.pie(fracs, labels=labels, shadow = True)
        plt.show()

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
