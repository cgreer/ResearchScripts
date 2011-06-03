import bioLibCG
import gZipEntropy
import cgNexusFlat
import cgOriginRNAFlat

def updateGScore(oFN, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['gScore', 'sequence'], [rn, tn])

        for oID in oNX.gScore:
                oNX.gScore[oID] = gZipEntropy.gZipEntropy(oNX.sequence[oID])


        oNX.save()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
