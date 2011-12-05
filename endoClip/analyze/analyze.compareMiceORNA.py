import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast
from compareData import compareTwoTcc
import cgOriginRNAFlat

def compareTwoMouse(aFN, bFN):

        aDC = cgNexusFlat.Nexus(aFN, cgOriginRNAFlat.OriginRNA)
        aDC.load(['tcc'])

        bDC = cgNexusFlat.Nexus(bFN, cgOriginRNAFlat.OriginRNA)
        bDC.load(['tcc'])
        
        aTccs = []
        bTccs = []

        for id in aDC.tcc:
                aTccs.append(aDC.tcc[id])

        for id in bDC.tcc:
                bTccs.append(bDC.tcc[id])

        
        print compareTwoTcc(aTccs, bTccs, amount = True)


def compareMouseHuman(mFN, hFNLift):

        aDC = cgNexusFlat.Nexus(mFN, cgOriginRNAFlat.OriginRNA)
        aDC.load(['tcc'])
        
        bTccs = []
        f = open(hFNLift, 'r')
        for line in f:
                ls = line.strip().split('\t')
                bTccs.append(ls[0])

        aTccs = []

        for id in aDC.tcc:
                aTccs.append(aDC.tcc[id])

        
        print compareTwoTcc(aTccs, bTccs, amount = True)


if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
