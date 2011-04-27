import bioLibCG
from siRnaPredict import getEntropy 


def filterEntropy(fN, outFN, minEntropy = 1.15):
        minEntropy = float(minEntropy)

        fOut = open(outFN, 'w')
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                seq = ls[1]
                ent = getEntropy(seq)

                if ent > minEntropy:
                        fOut.write(line)
                        

        

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(filterEntropy, sys.argv)


