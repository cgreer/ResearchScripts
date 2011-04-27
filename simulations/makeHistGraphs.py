import bioLibCG
import matplotlib.pyplot as plt
import siRnaPredict as si

def makeComplexityHist(fN):
        
        f = open(fN, 'r')

        histVals = []
        for line in f:
                ls = line.strip().split('\t')
                name = ls[3]
                seq = ls[4]

                if 'hsa' in name:
                        histVals.append(si.getEntropy(seq))

        plt.hist(histVals, 30)
        plt.show()


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(makeComplexityHist, sys.argv)


