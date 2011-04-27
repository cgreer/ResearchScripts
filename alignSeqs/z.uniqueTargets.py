import cgDB
import cgOriginRNA
import bioLibCG

def uniqueTargets(oDir):
        
        oDC = cgDB.dataController(oDir, cgOriginRNA.OriginRNA)
        id_oRNA = oDC.load()

        aID_numHit = {}
        uniqueTargets = []
        totalTargets = []
        for oRNA in id_oRNA.values():
                if not oRNA.passedFilter:
                        continue
                for aID in oRNA.filteredTargets:
                        numHits = aID_numHit.get(aID, 0)
                        aID_numHit[aID] = numHits + 1

                        if aID not in uniqueTargets:
                                uniqueTargets.append(aID)
                        
                        totalTargets.append(aID)

        for aID, numHit in aID_numHit.items():
                print aID, numHit

        print len(uniqueTargets)                
        print len(totalTargets)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(uniqueTargets, sys.argv)
