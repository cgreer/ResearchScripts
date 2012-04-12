import bioLibCG
from cgNexus import Nexus
import cgDL
from cgAutoCast import autocast
from cgAutoKeyWord import autokey
from bioLibJA import subit

def pickBestAlignment(fN, fFN):

    NX = Nexus(fN, fFN)
    NX.load(['sID', 'dID', 'score', 'best'])

    #find best id
    pair_score = {} #pair : score
    pair_highID = {}
    while NX.nextID():
        pair = '%s_%s' % (NX.sID, NX.dID)
        score = NX.score

        if score > pair_score.get(pair, 0.0):
            pair_score[pair] = score
            pair_highID[pair] = NX.id
   
    # update best id
    bestIDs = set(pair_highID.values())
    while NX.nextID():
        
        if NX.id in bestIDs:
            NX.best = True
    
    NX.save()

def updateIndexLists(qString, rString, qLen, rLen, qRange, rRange):
    '''used for getting the corresponding NT position on the reference
    strand from the query strand.  Used for peak center/mismatch center'''
    
    #update indexes
    currentQ = qRange[-1] # 3' end of query
    currentR = rRange[0] # 5' end of query
    qPositions = []
    rPositions = []
    for i, charMatch in enumerate(zip(qString, rString)):
       
        #query
        if '-' not in charMatch[0]:
            qPositions.append(currentQ)
            currentQ -= 1
        else:
            qPositions.append('-')

        #reference
        if '-' not in charMatch[1]:
            rPositions.append(currentR)
            currentR += 1
        else:
            rPositions.append('-')

    return (qPositions, rPositions)

def updateScores(fN, fFN):
    
    NX = Nexus(fN, fFN)
    NX.load(['numNormMatches', 'numGUs', 'numMismatches', 'numQGaps', 'numRGaps', 'numExtensionsQ','numExtensionsR', 'score'])
    
    while NX.nextID():
        NX.score = calculateAlignmentScore(NX.numNormMatches, NX.numGUs, NX.numMismatches, NX.numQGaps, NX.numRGaps, NX.numExtensionsQ, NX.numExtensionsR)

    NX.save()

def filterCenterProperties(fN, fFN):
    
    NX = Nexus(fN, fFN)
    NX.load(['query', 'reference', 'qStart', 'qEnd', 'rStart', 'rEnd', 'qLen', 'rLen', 'sigMask', 'centerPass', 'mismatchPass'])
   
    while NX.nextID():
        qRange = [NX.qStart, NX.qEnd]        
        rRange = [NX.rStart, NX.rEnd]        
        
        NX.mismatchPass = checkMismatchCenter(NX.query, NX.reference, NX.qLen, NX.rLen, qRange, rRange, NX.sigMask)
        NX.centerPass = checkPeakCenter(NX.query, NX.reference, NX.qLen, NX.rLen, qRange, rRange)

    NX.save()

def checkMismatchCenter(qString, rString, qLen, rLen, qRange, rRange, sigMask, middleRange=[9, 12]):
    '''given strings with gaps for query/reference, check 
    if mm will be in 9-12 nt of small RNA
    True if passed test
    NOTE: it seems like this takes care of the fact that the mask is reversed.  The alignments
    in the file have the oRNA reversed as well'''

    qPositions, rPositions = updateIndexLists(qString, rString, qLen, rLen, qRange, rRange)

    #find 9th and 12 position respective to alignment string reference
    indexOfFirst = qPositions.index(middleRange[1]) #12 first
    indexOfLast = qPositions.index(middleRange[0]) #now 9
    
    sigMask = list(sigMask)
    mmCheck = ('X' not in sigMask[indexOfFirst:indexOfLast + 1]) and ('G' not in sigMask[indexOfFirst:indexOfLast + 1])

    return mmCheck

def checkPeakCenter(qString, rString, qLen, rLen, qRange, rRange, middleRange=[9, 12]):
    '''given strings with gaps for query/reference, check 
    if peak will be in 9-12 nt of small RNA'''

    #init
    pPosition = rLen - 25

    #quick check if range even has a chance (peak not in range)
    if not pPosition in range(rRange[0], rRange[1]):
        return False
    
    ##aligned past edge of reference
    if rRange[0] == 1 or rRange[1] == rLen:
        return False

    qPositions, rPositions = updateIndexLists(qString, rString, qLen, rLen, qRange, rRange)

    #check position
    indexOfFirst = qPositions.index(middleRange[1]) #12 first
    firstRefNT = rPositions[indexOfFirst]

    indexOfLast = qPositions.index(middleRange[0]) #now 9
    secondRefNT = rPositions[indexOfLast]

    #print qPositions
    #print rPositions
    #print qString
    #print rString
    #print indexOfFirst, firstRefNT
    #print indexOfLast, secondRefNT
    #print pPosition
    #finally, check if peak is in center
    pCheck = (pPosition in rPositions[indexOfFirst:indexOfLast + 1])
    return pCheck

#TAG::alignment score, miRanda::
def calculateAlignmentScore(numM, numGU, numMM, numQGaps, numRGaps, numExtensionsQ, numExtensionsR):
    score = (5*numM) + (2*numGU) - (3*numMM) - (8 * (numQGaps + numRGaps)) - (2 * (numExtensionsQ + numExtensionsR))
    return score

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

