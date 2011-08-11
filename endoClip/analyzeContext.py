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
        oNX.load(['tccs', 'snrSS'], [rn, tn])

        histVals = []
        for oID in oNX.tccs:
                if not oNX.snrSS[oID] >= 2.0: continue
                histVals.append(len(oNX.tccs[oID]))

        plt.title('Number of Mapping Locations')
        plt.ylabel('Number of oRNA')
        plt.xlabel('Number of Places Mapped')                             
        plt.hist(histVals)
        plt.show()

def oRNAContextPie(oFN, rn = None, tn = None):
        '''REMEMBER!!! Have to do with grouped results...''' 
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['context', 'snrSS'], [rn, tn])

        context_count = {}

        for oID in oNX.context:
                
                if oNX.snrSS[oID] < 2.00: continue
                con = oNX.context[oID]
                print con
                context_count[con] = context_count.get(con, 0) + 1

        labels = sorted(context_count.keys())
        counts = [context_count[x] for x in labels]
        fracs = pieFractions(counts)

        #add numbers to labels
        labels = ['%s (%s)' % (x, context_count[x]) for x in labels]
        
        #plot
        plt.title('Context of oRNA (results > 2.00 SNR)')
        plt.pie(fracs, labels=labels, shadow = True)
        plt.show()

def oRNATypePie(oFN, rn = None, tn = None):        
        
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['snrSS', 'transcriptType', 'transcriptTypes'], [rn, tn])

        context_count = {}

        for oID in oNX.transcriptType:

                #con = oNX.transcriptType[oID]
                if oNX.snrSS[oID] < 2.00: continue
                cons = oNX.transcriptTypes[oID]
                for con in cons:
                        context_count[con] = context_count.get(con, 0) + 1
                #context_count[con] = context_count.get(con, 0) + 1

        labels = sorted(context_count.keys())
        counts = [context_count[x] for x in labels]
        fracs = pieFractions(counts)

        #add numbers to labels
        labels = ['%s (%s)' % (x, context_count[x]) for x in labels]
        
        #plot
        plt.pie(fracs, labels=labels, shadow = True)
        plt.show()

def targetContextPie(oFN, aFN, oContext = None, oType = None, rn = None, tn = None):
        if oContext == 'None': oContext = None
        if oType == 'None': oType = None
        
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['snrSS', 'context', 'transcriptType', 'filteredTargets'], [rn, tn])

        aNX = cgNexusFlat.Nexus(aFN, cgAlignmentFlat.cgAlignment)
        aNX.load(['context'], [rn, tn])

        context_count = {}

        for oID in oNX.context:

                #filter oID types here
                con = oNX.context[oID]
                typ = oNX.transcriptType[oID]
                
                if oContext:
                        if oContext != con: continue

                if oType:
                        if oType != typ: continue

                
                if oNX.snrSS[oID] < 2.00: continue
                #gather targets' context info if oRNA is okay
                for aID in oNX.filteredTargets[oID]:
                        aCon = aNX.context[aID]
                        if aCon == 'C_3UTR': print oID, aID
                        context_count[aCon] = context_count.get(aCon, 0) + 1

        labels = sorted(context_count.keys())
        counts = [context_count[x] for x in labels]
        fracs = pieFractions(counts)

        #add numbers to labels
        labels = ['%s (%s)' % (x, context_count[x]) for x in labels]
        
        #plot
        plt.title('Context of Targets (results > 2.00 SNR)\n oContext: %s, oType: %s' % (oContext, oType))
        plt.pie(fracs, labels=labels, shadow = True)
        plt.show()

def targetTypePie(oFN, aFN, oContext = None, oType = None, rn = None, tn = None):
        if oContext == 'None': oContext = None
        if oType == 'None': oType = None


        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['filteredTargets','transcriptType', 'context'], [rn, tn])

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
        plt.title('Type of Targets (results > 2.00 SNR)\n oContext: %s, oType: %s' % (oContext, oType))
        plt.pie(fracs, labels=labels, shadow = True)
        plt.show()

def targetContextPercentageVsExpression(oFN, aFN, oContext = None, oType = None, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['context', 'transcriptType', 'filteredTargets'])

        aNX = cgNexusFlat.Nexus(aFN, cgAlignmentFlat.cgAlignment)
        aNX.load(['context', 'tELevel'])

        context_level_count = {}

        for lev in range(50, 500, 50):
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
                                eLevel = aNX.tELevel[aID]
                                if eLevel < lev: continue
                                aCon = aNX.context[aID]
                                context_level_count[aCon][lev] = context_level_count.setdefault(aCon, {}).get(lev, 0) + 1


        #fracs = pieFractions(counts)
        plots_labels = [[], []]
        for con in context_level_count:
                x = []
                y = []
                sortedLevs = sorted(context_level_count[con].keys())
                for lev in sortedLevs:
                        x.append(lev)
                        y.append(context_level_count[con][lev])
                plots_labels[0].append(plt.plot(x,y,))
                plots_labels[1].append(con)
                
                                        
        
        #plot
        plt.legend(plots_labels[0], plots_labels[1])
        plt.title('oRNA Targets\' Context Proportion Stability w/ Expression Increase')
        plt.xlabel('Degradome Expression Cutoff')
        plt.ylabel('Number of oRNA Targets')
        plt.show()

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
