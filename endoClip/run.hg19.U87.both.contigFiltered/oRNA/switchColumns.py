import bioLibCG
from cgNexus import Nexus
import cgDL
from cgAutoCast import autocast
from cgAutoKeyWord import autokey
from bioLibJA import subit

def switchColumns(fN, outFN):

    fOut = open(outFN, 'w')
    f = open(fN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        ls[2], ls[3] = ls[3], ls[2]
        fOut.write('\t'.join(ls) + '\n')

    f.close()
    fOut.close()

    

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

