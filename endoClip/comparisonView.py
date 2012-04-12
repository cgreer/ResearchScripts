import bioLibCG
from cgNexus import Nexus
import cgDL
from cgAutoCast import autocast
from cgAutoKeyWord import autokey
from bioLibJA import subit

def outputExpandedAlignment(oFN, oFF, dFN, dFF, aFN, aFF, outFN):

    oNX = Nexus(oFN, oFF)
    oNX.load(['tcc', 'eLevel', 'filteredTargets', 'avgNumSS', 'snr', 'numUFBS'])

    dNX = Nexus(dFN, dFF)
    dNX.load(['tcc', 'eLevel', 'pValBin'])

    aNX = Nexus(aFN, aFF)
    aNX.load(['sID', 'dID', 'query', 'reference', 'sigMask', 'adjustedNumMismatches'])

    with open(outFN, 'w') as f:
        while oNX.nextID():
            for aID in oNX.filteredTargets:
                aNX.id = aID
                dNX.id = aNX.dID
                outS = [oNX.id, dNX.id, aNX.id, oNX.tcc, dNX.tcc, oNX.eLevel, dNX.eLevel, oNX.numUFBS, 

    

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

