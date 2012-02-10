import bioLibCG
from cgAutoCast import autocast


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

def removeContigs(fN, outFN):
    
    readPass = False
    currentRead = [] 
    outF = open(outFN, 'w')
    f = open(fN, 'r')
    for i, line in enumerate(f):

        #bioLibCG.printSeldomly('processed %s lines' % i, 20)

        if (i - 1) % 4 == 0:
            #check contigs
            seq = line.strip()
            m = getLongestRepeat(seq, 4, 1) 
            d = getLongestRepeat(seq, 2, 2) 
            t = getLongestRepeat(seq, 2, 3) 
            q = getLongestRepeat(seq, 2, 4) 
            
            #print m, d, t, q
            #print seq
            readPass = True
            if (m > 5) or (d > 3) or (t > 2) or (q > 2):
                readPass = False
            
            currentRead.append(line)
        elif (i -3) % 4 == 0:
            #save read
            if readPass:
                currentRead.append(line)
                outF.writelines(currentRead)
            
            currentRead = []
            
        else:
            currentRead.append(line)
    
    f.close()
    outF.close()


if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
