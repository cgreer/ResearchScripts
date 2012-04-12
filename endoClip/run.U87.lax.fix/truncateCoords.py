import bioLibCG
from cgNexus import Nexus
import cgDL
from cgAutoCast import autocast
from cgAutoKeyWord import autokey
from bioLibJA import subit
from cgLog import Logger 
logger = Logger(0)

@autocast
def truncate(fN, spread = 3):
   
    f = open(fN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        chrom, strand, start, end  = ls[0].split(":")
        for i in range(1,spread + 1):
            print '%s:%s:%s' % (chrom, strand, int(start) + i)
            print '%s:%s:%s' % (chrom, strand, int(start) - i)

    f.close()
    

if __name__ == "__main__":
    import sys
    assert sys.argv[1] in globals(), "Need name of fxn to run from command line!"
    fxnToRun = globals()[sys.argv[1]] 
    fxnToRun(*sys.argv[2:])

