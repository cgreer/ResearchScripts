import bioLibCG
from cgNexus import Nexus
import cgDL
from cgAutoCast import autocast
from cgAutoKeyWord import autokey
from bioLibJA import subit
import math

@autocast
def collectData(fN, fFN, outFN, logHeight = False):
    '''x(pHeight) v. y(-log10(pval)).  0.0 is 1.0e-100'''

    NX = Nexus(fN, fFN)
    NX.load(['eLevel', 'pValBin'])

    f = open(outFN, 'w')
    while NX.nextID():

        x = math.log(NX.eLevel, 10) if logHeight else NX.eLevel
        pVal = NX.pValBin
        if pVal < 0: continue
        pVal = pVal if (pVal != 0.0) else float("1.0e-100") 
                
        try:
            y = -math.log(pVal, 10)
        except ValueError:
            print x, pVal
            return

        f.write('%s\t%s\n' % (x,y))
    f.close()
    

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

