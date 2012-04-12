import bioLibCG
from cgNexus import Nexus
import cgDL
from cgAutoCast import autocast
from cgAutoKeyWord import autokey
from bioLibJA import subit

def getMicroHistExpression(microFN, fqFile, outFN):

    microSeqs = cgDL.listFromColumns(microFN, [0], ['string'])
    microSeq_count = dict( (seq, 0) for seq in microSeqs)
        
    f = open(fqFile, 'r')
    for line in f:
        possibleSeq = line.strip()
        if possibleSeq in microSeq_count:
            microSeq_count[possibleSeq] += 1
    f.close()

    with open(outFN, 'w') as f:
        for seq, count in microSeq_count.iteritems():
            f.write('%s\t%s\n' % (seq, count))
    

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

