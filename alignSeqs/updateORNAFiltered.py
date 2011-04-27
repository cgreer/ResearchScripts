import cgNexusFlat
import cgOriginRNAFlat
import bioLibCG

def updateFiltered(oFN):
        
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['filteredTargets', 'endContigLength', 'totalContigLength', 'entropy', 'sequenceDuplicate', 'passedFilter'])

        for oID in oNX.passedFilter:

                oNX.passedFilter[oID] = False

                if len(oNX.filteredTargets[oID]) == 0:
                        continue

                if oNX.endContigLength[oID] > 6:
                        continue

                if oNX.totalContigLength[oID] > 6:
                        continue

                if oNX.entropy[oID] < 1.2:
                        continue

                if oNX.sequenceDuplicate[oID]:
                        continue

                #if it passed, update
                oNX.passedFilter[oID] = True




        oNX.save()
        

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
