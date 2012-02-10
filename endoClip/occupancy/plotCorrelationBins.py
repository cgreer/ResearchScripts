import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast
import matplotlib.pyplot as plt

def plotCorrelationBins(fN):
    '''hist only works when there is one column...
    multi columns works only for the same file...'''

    f = open(fN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        Xs = [float(x) for x in ls[2].split(',')]
        Ys = [float(x) for x in ls[3].split(',')]
    f.close()

    pairs = zip(Xs, Ys)
    newPairs = [pair for pair in pairs if (pair[0] != 0 and pair[1] != 0)]
    newXs = [pair[0] for pair in newPairs]
    newYs = [pair[1] for pair in newPairs]

    
        
    plt.scatter(newXs, newYs)
    plt.show()
    

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

