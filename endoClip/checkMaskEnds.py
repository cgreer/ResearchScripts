import bioLibCG
from cgNexus import Nexus
import cgDL
from cgAutoCast import autocast
from cgAutoKeyWord import autokey
from bioLibJA import subit

def checkMaskEnds(maskPerLineFN):

    masks = cgDL.listFromColumns(maskPerLineFN, [0], ['string'])

    index_numMM = dict((i,0) for i in range(10))

    for mask in masks:
        mask = mask[::-1]
        for i, char in enumerate(mask):
            if i == 10: break
            if char == 'X':
                index_numMM[i] += 1


    for i, num in index_numMM.items():
        print i, num

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

