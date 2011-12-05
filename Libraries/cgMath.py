import bioLibCG
from cgAutoCast import autocast

@autocast
def linspace(bottomRange, topRange, numPieces):
    '''Return evenly spaced numbers over a specified interval
    Same as numpy version --> needed to run pypy without numpy'''

    rangeLength = (float(topRange) - float(bottomRange))
    spacer = rangeLength / (numPieces - 1)
    linSpace = [float(bottomRange)]
    ext = [spacer * x for x in range(1, numPieces - 1)]
    linSpace.extend(ext)
    linSpace.append(float(topRange))

    return linSpace
    

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
