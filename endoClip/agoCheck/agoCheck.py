import bioLibCG
from cgNexus import Nexus
import cgDL
from cgAutoCast import autocast
from cgAutoKeyWord import autokey
from bioLibJA import subit

def nextFilePacket(fHandler, linesPerPacket):
    '''useful for fa and fastq etc'''

    nextPacket = [fHandler.readline() for x in range(linesPerPacket)]

    emptyCheck = [x == '' for x in nextPacket]
    if all(emptyCheck):
        return None
    else:
        return nextPacket

@autocast
def check_ORNA_in_ago(oFN, oFF, agoFN, clippingAmount = 1):
   
    NX = Nexus(oFN, oFF)
    NX.load(['sequence', 'geneNames'])
    
    #make truncated sequences
    id_sequence = NX.createMap('id', 'sequence')
    if clippingAmount > 0:
        id_sequence = dict( (i, j[clippingAmount:-clippingAmount]) for i,j in id_sequence.items())
   
    #get fastq sequences
    agoF = open(agoFN, 'r')
    agoSeqs = []
    while True:
        fPacket = nextFilePacket(agoF, 4)
        if not fPacket: break
        agoSeqs.append(fPacket[1])
    agoF.close()

    #count for each oRNA
    id_count = {}
    for id, seq in id_sequence.items():
        for agoSeq in agoSeqs:
            if seq in agoSeq:
                id_count[id] = id_count.get(id, 0) + 1

    #out
    totalCount = 0
    for id, count in id_count.items():
        NX.id = id
        print '%s\t%s\t%s' % (id, count, NX.geneNames)
        totalCount += count

    print totalCount
    

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
