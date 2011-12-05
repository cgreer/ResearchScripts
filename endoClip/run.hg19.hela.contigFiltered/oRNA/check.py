import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast
from cgOriginRNAFlat import OriginRNA

def numTargets(fN):
    
    NX = cgNexusFlat.Nexus(fN, OriginRNA)
    NX.load(['filteredTargets'])

    tSet = set()
    ts = []
    for id in NX.ids:
        for tID in NX.filteredTargets[id]:
            tSet.add(tID)
            ts.append(tID)

    print len(tSet)
    print len(ts)
            


if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
