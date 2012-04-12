import bioLibCG
import cgNexusFlat
import cgDL
from cgAutoCast import autocast
import string
from cgNexus import Nexus

@autocast
def updateAdjustedMismatches(fN, fF, guValue = .5, otherValue = 1.0):
    
    NX = Nexus(fN, fF)
    NX.load(['sigMask', 'adjustedNumMismatches'])

    while NX.nextID():

        mask = NX.sigMask
        numGU = mask.count('G')
        numGapAndMM = mask.count('X')
        
        NX.adjustedNumMismatches = (numGU * guValue) + (numGapAndMM * otherValue)

    NX.save()

@autocast
def updateAdjustedMismatchesFlat(fN, fOut, guValue = .5, otherValue = 1.0):
    
    f = open(fN, 'r')
    fOut = open(fOut, 'w')
    for line in f:
        ls = line.strip().split('\t')
        
        mask = ls[18] 
        numGU = mask.count('G')
        numGapAndMM = mask.count('X')
        
        adjustedNumMismatches = (numGU * guValue) + (numGapAndMM * otherValue)
        ls.append(str(adjustedNumMismatches))
        fOut.write('\t'.join(ls) + '\n') 
    f.close()
    fOut.close()

#TAG::lowercase,count from end::
def countLowerEnd(theString, fromEnd = False):

    if fromEnd:
        theString = theString[::-1]

    numLower = 0
    for c in theString:
        if c in string.lowercase:
            numLower += 1
        else:
            return numLower

#TAG::collapse,repeat,string::
def collapseRuns(theString):

    prevChar = None
    newString = []
    for c in theString:
        if c == prevChar:
            continue
        else:
            newString.append(c)
            prevChar = c

    return ''.join(newString)

#TAG::read multiple file lines,read file,file::
def lineFileParser(inFile, linesPerEntry = 4):
    '''retrieve X lines at a time'''

    f = open(inFile, 'r')

    allLines = []
    while True:
        newEntry = [f.readline() for i in range(linesPerEntry)]
        
        for l in newEntry:
            if l == '':
                return allLines

        allLines.append(newEntry)

def parseRandaOutput(inFile, smallFN, degFN, oFN):
    '''miRanda Output will be grepped to make it every alignment is 10lines'''

    allAlignments = lineFileParser(inFile, 10)
    smallID_length = getIDLength(smallFN)
    degID_length = getIDLength(degFN)

    outF = open(oFN, 'w')
    for i, alignment in enumerate(allAlignments):
        oldScore, qRange, rRange, query, matchInfo, reference, sID, dID = parseAlignmentBasic(alignment)
        complexResult = parseAlignmentComplex(query, reference)
        numMM, numM, numGU, numQGaps, numRGaps, significantTargetMask, numExtensionsQ, numExtensionsR = complexResult
        
        #check for alignments with N in them
        if all([x == 0 for x in complexResult]): continue

        sLen = smallID_length[int(sID)]
        dLen = degID_length[int(dID)]

        pString = [i, sID, dID, qRange[0], qRange[1], rRange[0], rRange[1], sLen, dLen, query, reference, numM, numMM, numGU, numQGaps, numRGaps, numExtensionsQ, numExtensionsR, significantTargetMask]
        pString = '\t'.join([str(x) for x in pString])
        outF.write(pString + '\n')
    outF.close()
    
#TAG::aligning,parse alignment,miRanda::
def parseAlignmentBasic(alignment):

    #parse raw data
    info, n1, query, matchInfo, reference, n3, n4, n5, idInfo, n6 = alignment

    oldScore = float(info.split()[2])

    qRange = (int(info.split()[3][2:]), int(info.split()[5]))
    rRange = (int(info.split()[6][2:]), int(info.split()[8]))

    query = query.split()[2]
    matchInfo = matchInfo.strip()
    reference = reference.split()[2]

    sID = idInfo.split()[0][1:]
    dID = idInfo.split()[1]

    #calculate read qRange and rRange
    '''qrange is weird, the lower number is correct and the higher one is to high by one (2 to 14) should be (2 to 13)
    In addition, the qRange has to have the lower case letters added to it whereas the rRange already includes it'''
    qRange = qRange[0] - countLowerEnd(query, fromEnd = True), qRange[1] + countLowerEnd(query) - 1 #should always be (1,N)  
    
    output = [oldScore, qRange, rRange, query, matchInfo, reference, sID, dID]
    if any([x == '' for x in output]):
        print output
        raise NameError("missing some parsing info")
    return output 

def parseAlignmentComplex(query, reference):
    '''Get num missmatches and other parsed data from alignment'''

    allLetters = ['A', 'T', 'C', 'G']
    queryGapPairs = ['-%s' % x for x in allLetters]
    referenceGapPairs = ['%s-' % x for x in allLetters]

    #TAG::genomic letter combinations,combinations,letters::
    guPairs = ['GT', 'TG']
    compPairs = ['AT', 'TA', 'CG', 'GC']
    misPairs = ['TT', 'TC', 'AG', 'AC', 'AA', 'CC', 'CT', 'CA', 'GG', 'GA'] 
    
    if len(reference) != len(query):
        raise NameError("ALIGNMENTS ARE DIFFERENT SIZE!")

    #calculate # gaps// ALLOW TARGET GAPS???
    collQ = collapseRuns(query)
    collR = collapseRuns(reference)
    numQGaps = collQ.count('-')
    numRGaps = collR.count('-')
    numExtensionsQ = query.count('-') - numQGaps
    numExtensionsR = reference.count('-') - numRGaps

    #calc match/mismatch (rev to get from small 5-->3)
    query = query.upper()[::-1]
    reference = reference.upper()[::-1]
    matchPairs = ['%s%s' % (x,y) for x,y in zip(query, reference)]

    significantTargetMask = [] #Mask is entire small string masked
    numM = numMM = numGU = 0
    gapShift = 0
    for i,pair in enumerate(matchPairs):
        if pair in queryGapPairs:
            gapShift -= 1
            numMM += 1
            significantTargetMask.append('X')
        elif pair in guPairs:
            numGU += 1
            significantTargetMask.append('G')
        elif pair in compPairs:
            numM += 1
            significantTargetMask.append('N')
        elif pair in misPairs:
            numMM += 1
            significantTargetMask.append('X')
        elif pair in referenceGapPairs:
            numMM += 1
            significantTargetMask.append('X')
        elif 'N' in pair:
            return [0,0,0,0,0,0,0,0] #dont take N alignments
        else:
            print query
            print reference
            print pair
            raise NameError("COMBINATION NOT ACCOUNTED FOR!!!")

    significantTargetMask = ''.join(significantTargetMask)
    significantTargetMask = significantTargetMask[::-1] # did re-aligning reverse from miRanda...for sanity
    return [numMM, numM, numGU, numQGaps, numRGaps, significantTargetMask, numExtensionsQ, numExtensionsR]

def getIDLength(fN):

    allSeqs = lineFileParser(fN, 3)
    id_length = {}
    for small in allSeqs:
        id, seq, blank = small
        id = int(id.strip().split('>')[-1])
        id_length[id] = len(seq.strip())

    return id_length

def parseRandaInclusiveCheck(inFile, smallFN):
    '''miRanda Output will be grepped to make it every alignment is 10lines'''

    id_length = getIDLength(smallFN)

    allAlignments = lineFileParser(inFile, 10)
    for i, alignment in enumerate(allAlignments):
        checkInclusiveSmallLength(alignment, id_length)

def checkInclusiveSmallLength(alignment, id_length):
    '''hacked script to check if miRanda always shows the full length small
    RNA at the QUERY part of the alignment...needed for alignment calculations'''

    #parse raw data
    info, n1, query, matchInfo, reference, n3, n4, n5, idInfo, n6 = alignment

    oldScore = float(info.split()[2])

    qRange = (int(info.split()[3][2:]), int(info.split()[5]))
    rRange = (int(info.split()[6][2:]), int(info.split()[8]))

    query = query.split()[2]
    matchInfo = matchInfo.strip()
    reference = reference.split()[2]

    sID = idInfo.split()[0][1:]
    dID = idInfo.split()[1]

    query = list(query)
    dashCount = query.count('-')
    if len(query) - dashCount != id_length[int(sID)]:
        print sID, dID
        raise NameError("MISMATCH")

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

