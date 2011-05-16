import bioLibCG
import cgNexusFlat
import cgOriginRNAFlat
import cgAlignmentFlat

def updateSignificantSequences(oFN, aFN, rn = None, tn = None):

        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['filteredTargets', 'numSignificantSequences'], [rn, tn])
        
        aNX = cgNexusFlat.Nexus(aFN, cgAlignmentFlat.cgAlignment)
        aNX.load(['mismatchPositions'])

        for oID in oNX.filteredTargets:

                mismatchSet = set()
                
                for aID in oNX.filteredTargets[oID]:

                        mismatchSet.add(tuple(aNX.mismatchPositions[aID]))

                oNX.numSignificantSequences[oID] = len(mismatchSet)

        oNX.save()

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
