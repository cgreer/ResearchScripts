import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast

def newCoords(fN):
        
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                chrom, strand, start, end = bioLibCG.tccSplit(ls[0])
                print '%s:%s-%s' % (chrom, start, end)

def toTcc(fN, strand):

        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                
                chrom = ls[0].split(':')[0]
                start = ls[0].split(':')[1].split('-')[0]
                end = ls[0].split(':')[1].split('-')[1]

                newTcc = bioLibCG.makeTcc(chrom, strand, start, end)
                print newTcc
                

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
