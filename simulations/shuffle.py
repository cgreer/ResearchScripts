import random
import bioLibCG
from cgNexus import Nexus
import cgDL
from cgAutoCast import autocast
from cgAutoKeyWord import autokey
from endoShuffle import shuffleMaskedSeq

def kMask(seq, mask, minTimes, kmerLength):

    if not mask:
        mask = [False for x in seq]
    #get all di seqs, uniquify
    kSeqs = bioLibCG.returnFrames(seq, kmerLength)
    kSeqs = set(kSeqs)


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
                            maskStart = i - (kmerLength * sLen)
                            for i in range(maskStart, i): mask[i] = True
                            sLen = 0
                        else:
                            sLen = 0
                except IndexError:
                    #check for masking
                    if sLen > minTimes:
                        maskStart = i - (kmerLength * sLen)
                        for i in range(maskStart, i): mask[i] = True
                        sLen = 0
    
    return mask

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

def maskAll(fN):

    allSeqs = cgDL.listFromColumns(fN, [2], ['string'])
    allMasks = [multiMask(x, list(), True) for x in allSeqs]
    

def multiMask(seq, mask = [], toPrint = False):
     
        toPrint = True if toPrint == 'True' else False
        mask = mask if mask != '.' else []
        
        mask = kMask(seq, mask, 5, 1)
        #print '   mask 1', mask
        mask = kMask(seq, mask, 3, 2)
        #print '   mask 2', mask
        mask = kMask(seq, mask, 2, 3)
        #print '   mask 3', mask
        mask = kMask(seq, mask, 2, 4)
        #print '   mask 4', mask
        mask = kMask(seq, mask, 2, 5)
        #print '   mask 5', mask

        if toPrint:
            print ''.join(['X' if x else '_' for x in mask])
        return mask

class Sequence:
        '''NOTE THAT THE SHUFFLE(not shufflemer) FXN WILL DISRUPT THE MER POSITIONS...'''

        def __init__(self, sequence, mer = None):
                self.sequence = sequence
                self.mer = mer
                self.merStarts = []
                self.merPositions = []
                self.merUpdateFlag = False
                self.findMerPositions()
     
        def getNucleotideFrequency(self):
                l = len(self.sequence)
                freqDict = {}
                for nt in ['A', 'C', 'T', 'G']:
                        num = self.sequence.count(nt)
                        freqDict[nt] = float(num)/l

                return freqDict

        def findMerPositions(self):
                l = len(self.sequence)
                lmer = len(self.mer)

                for i in range(0, l):
                        if self.sequence[i:i + lmer] == self.mer:
                                self.merPositions.extend([i,i + 1])
                                self.merStarts.append(i)
                
                self.merUpdateFlag = True

        def shuffleNucleotide(self, intactMer = False):

                seqL = len(self.sequence)
                merLength = len(self.mer)

                          
                #first choose a starting mer/nt to swap
                merStart = random.randint(0, seqL - 1)                
                swapPosition = random.randint(0, seqL - 1)# pick so that the end of mer is not beyond end of sequence

                #check to see if the swap position is off limits
                if intactMer:
                        if (swapPosition) in self.merPositions or merStart in self.merPositions:
                                return False

                #if it passed, then go ahead and make the swap + change the mer positions
                self.sequence = list(self.sequence)
                self.sequence[merStart], self.sequence[swapPosition] = self.sequence[swapPosition], self.sequence[merStart]
                self.sequence = ''.join(self.sequence)
                
                return True

        def shuffleMer(self):
                '''Note: This doesn't work on overlapping mers such as AAA in AAAA.  The mer will be destroyed...'''

                seqL = len(self.sequence)
                merLength = len(self.mer)

                          
                #first choose a starting mer/nt to swap
                merStart = random.choice(self.merStarts)                
                swapPosition = random.randint(0, seqL - 1 - (merLength - 1))# pick so that the end of mer is not beyond end of sequence

                #check to see if the swap position is off limits
                for i in range(0, merLength):
                        if (swapPosition + i) in self.merPositions:
                                return False

                #if it passed, then go ahead and make the swap + change the mer positions
                for i in range(0, merLength):
                        self.sequence = list(self.sequence)
                        self.sequence[merStart + i], self.sequence[swapPosition + i] = self.sequence[swapPosition + i], self.sequence[merStart + i]
                        self.sequence = ''.join(self.sequence)
                        
                        if i == 0:
                                self.merStarts.remove(merStart)
                                self.merStarts.append(swapPosition)
                        
                        self.merPositions.remove(merStart + i)
                        self.merPositions.append(swapPosition + i)

                return True

def getTotalContigLength(seq):

        highestLength = 1
        cLength = 1
        letters = list(seq)
        for i,letter in enumerate(letters):
                if i == 0: continue

                if letters[i] == letters[i-1]:
                        cLength += 1
                        if cLength > highestLength:
                                highestLength = cLength
                else:
                        cLength = 1
        
        return highestLength

def seqShuffle(seq, mask, timesShuffled = 1):

    seqL = len(seq)
    numFalseShuffles = 0
    for i in xrange(0, timesShuffled):
            
        #first choose a starting mer/nt to swap
        ntStart = random.randint(0, seqL - 1)                
        ntEnd = random.randint(0, seqL - 1)# pick so that the end of mer is not beyond end of sequence

        #check to see if the swap position is off limits
        if mask[ntStart] or mask[ntEnd]:
            numFalseShuffles +=1
            continue


        #if it passed, then go ahead and make the swap + change the mer positions
        seq = list(seq)
        seq[ntStart], seq[ntEnd] = seq[ntEnd], seq[ntStart]
        seq = ''.join(seq)
    
    
    if numFalseShuffles/float(timesShuffled) > .90: print '..low efficiency shuffle:', seq, numFalseShuffles, timesShuffled
    return seq

def passContigFilter(seq):
    
    m = getLongestRepeat(seq, 4, 1) 
    d = getLongestRepeat(seq, 2, 2) 
    t = getLongestRepeat(seq, 2, 3) 
    q = getLongestRepeat(seq, 2, 4) 

    if (m > 5) or (d > 3) or (t > 2) or (q > 2):
        return False
    else:
        return True

def maskCGs(seq, mask = []):
    '''Think it's faster than regex even though a little crude'''

    if not mask:
        mask = [False for x in range(len(seq))]

    seqLength = len(seq)

    #first frame
    for i, let in enumerate(seq):
        if i == seqLength - 1: continue

        if seq[i:i+2] == 'CG':
            mask[i] = True
            mask[i + 1] = True
    
    
    #2nd frame
    for i, let in enumerate(seq[1:]):
        if i == seqLength - 2: continue

        if seq[i:i+2] == 'CG':
            mask[i + 1] = True
            mask[i + 2] = True

    return mask

@autocast
def shuffleSeq2(fN, outFN, fasta = True):
        
    f = open(fN, 'r')
    fOut = open(outFN, 'w')
    
    for line in f:
        ls = line.strip().split('\t')
        oID = ls[0]
        seq = ls[1]
        
        #have to convert 
        seq = seq.replace("T", "U")
       
        #mask the sequence for repeats, if everything is masked than mark bad egg
        mask = multiMask(seq, [])
        if all(mask):
            print oID, "is not shufflin'....everyday"
            continue
        
        #try to shuffle up to 100 times
        for i in range(100):
            shuffledSeq = shuffleMaskedSeq(seq, mask)
            if passContigFilter(shuffledSeq):

                #revert 
                shuffledSeq = shuffledSeq.replace("U", "T")

                #output correct format
                if fasta:
                    fOut.write('>%s\n%s\n\n' % (oID, shuffledSeq))
                else:
                    fOut.write('%s\t%s\n' % (oID, shuffledSeq))

                #no need to reshuffle
                break

            if i == 99:
                print oID, "is not shufflin'....everyday"

    f.close()
    fOut.close()

def shuffleSeq(fN, outFN):

        f = open(fN, 'r')
        fOut = open(outFN, 'w')
        
        for line in f:
                ls = line.strip().split('\t')
                sequence = ls[1]
                oID = ls[0]
                seq = Sequence(sequence, 'CG')
               
                while True:
                        
                        #shuffle mers
                        if 'CG' in sequence:
                                j = 0
                                while j < 50:
                                        if seq.shuffleMer() == False:
                                                continue
                                        else:
                                                j += 1
                        
                        #shuffle nts
                        j = 0
                        while j < 100:
                                if seq.shuffleNucleotide(True) == False:
                                        continue
                                else:
                                        j += 1
                        
                        if getTotalContigLength(seq.sequence) < 7:
                                fOut.write('%s\t%s\n' % (oID, seq.sequence))
                                break
                        else:
                                print 'did not pass contig filter'

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
