import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast

def consolidatePeaksByTcc(dFN, outFN):

        fOut = open(outFN, 'w')
        f = open(dFN, 'r')

        lastLine = ''
        lastTcc = 'chr100:1:100:1000'
        for line in f:
                ls = line.strip().split('\t')
                tcc = ls[1]

                
                
                if bioLibCG.tccOverlap(tcc, lastTcc):
                        #skip line...
                        lastLine = line
                        lastTcc = tcc
                else:
                        fOut.write(lastLine)
                        lastLine = line
                        lastTcc = tcc
                
        fOut.close()
        f.close()


if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
