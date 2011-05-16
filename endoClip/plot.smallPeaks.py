import bioLibCG
import cgNexusFlat
import cgOriginRNAFlat

def smallPeaks(fN, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(fN, className)
        DC.load(['tcc'], [rn, tn])

        


        oNX.save()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
