import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast
import microSeq


@autocast
def getLongestRepeat(seq, minTimes, kmerLength):
    
    #get all di seqs, uniquify
    kSeqs = bioLibCG.returnFrames(seq, kmerLength)
    kSeqs = set(kSeqs)
    seqLength = len(seq)

    highestSLen = 0
    for kmer in kSeqs:
        for slide in range(0, kmerLength):
            sLen = 0
            adjustedRange = seqLength - (kmerLength - 1) #takes into account excluding right side numbers
            scanRange = range(slide, adjustedRange, kmerLength)
            for i in scanRange:
                #print seq[:i] + ' [' + seq[i:i + kmerLength] + ']', highestSLen, sLen, i, i + kmerLength
                if (seq[i:i + kmerLength] == kmer):
                    sLen += 1
                    
                    #take care of last nt
                    if i == scanRange[-1]:
                        if sLen >= minTimes:
                            if sLen > highestSLen: highestSLen = sLen
                            sLen = 0

                else:
                    #if stretch is long enough, mask
                    if sLen >= minTimes:
                        if sLen > highestSLen: highestSLen = sLen
                    sLen = 0                            
    
    return highestSLen

def updateTotalContig(oFN, rn = None, tn = None):

    oNX = cgNexusFlat.Nexus(oFN, microSeq.MicroSeq)
    oNX.load(['seq', 'longestMono', 'longestDi', 'longestTri', 'longestQuad'], [rn, tn])

    for oID in oNX.ids:
            seq = oNX.seq[oID] 
            oNX.longestMono[oID] = getLongestRepeat(seq, 1, 1)
            oNX.longestDi[oID] = getLongestRepeat(seq, 1, 2)
            oNX.longestTri[oID] = getLongestRepeat(seq, 1, 3)
            oNX.longestQuad[oID] = getLongestRepeat(seq, 1, 4)
            
            
    oNX.save()

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
