import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast

def name(fN):
    
    f = open(fN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        
        seq = ls[3]
        numTargets = len(ls[9].split(','))
        print '%s\t%s' % (seq, numTargets)

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
