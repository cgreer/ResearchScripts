import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast
import microSeq

@autocast
def getLongestRepeat(seq, minTimes, kmerLength):
    
    #get all di seqs, uniquify
    kSeqs = bioLibCG.returnFrames(seq, kmerLength)
    kSeqs = set(kSeqs)

    highestSLen = 0
    for kmer in kSeqs:
        #slides
        for slide in range(0, kmerLength):

            sLen = 0
            for i in range(slide, len(seq), kmerLength):

                try:
                    if seq[i:i + kmerLength] == kmer:
                        sLen += 1
                    else:
                        #if stretch is long enough, mask
                        if sLen >= minTimes:
                            if sLen > highestSLen: highestSLen = sLen
                            sLen = 0
                        else:
                            sLen = 0
                except IndexError:
                    #check for masking
                    if sLen > minTimes:
                        if sLen > highestSLen: highestSLen = sLen
                        sLen = 0
    print highestSLen                        
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
