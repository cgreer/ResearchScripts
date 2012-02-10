import bioLibCG
import cgNexusFlat
import cgDL
from cgAutoCast import autocast
import string

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

def pickBestAlignment(inFile, outFN):

    pair_scoreLine = {} #pair : score, line
    f = open(inFile, 'r')
    for line in f:
        ls = line.strip().split('\t')
        pair = '%s_%s' % (ls[1], ls[2])
        try:
            score = float(ls[9])
        except:
            print line
            raise NameError("BLAH")

        if score > pair_scoreLine.get(pair, [0,0])[0]:
            pair_scoreLine[pair] = [score, line]
    f.close()
   
    f = open(outFN, 'w')
    for pair, [score, line] in pair_scoreLine.iteritems():
        f.write(line)
    f.close()
    
def parseRandaOutput(inFile, oFN):
    '''miRanda Output will be grepped to make it every alignment is 10lines'''

    allAlignments = lineFileParser(inFile, 10)

    outF = open(oFN, 'w')
    for i, alignment in enumerate(allAlignments):
        oldScore, qRange, rRange, query, matchInfo, reference, sID, dID = parseOneAlignment(alignment)
        aScore, numMM, numM, numGU, numGaps, smallMMPos = calculateAlignmentScore(query, reference)
        if aScore == 0: continue #takes care of alignments with N in them
        pString = [i, sID, dID, qRange[0], qRange[1], rRange[0], rRange[1], numMM, ','.join([str(x) for x in smallMMPos]), aScore, numGaps]
        pString = '\t'.join([str(x) for x in pString])
        outF.write(pString + '\n')
    outF.close()
    
#TAG::print by line::
def printBL(pInfo):
    for l in pInfo:
        print l.strip()

#TAG::aligning,parse alignment,miRanda::
def parseOneAlignment(alignment):

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
    qRange = qRange[0] - countLowerEnd(query, fromEnd = True), qRange[1] + countLowerEnd(query)  
    rRange = rRange[0] - countLowerEnd(reference), rRange[1] + countLowerEnd(reference, fromEnd = True)  
    
    output = [oldScore, qRange, rRange, query, matchInfo, reference, sID, dID]
    if any([x == '' for x in output]):
        print output
        raise NameError("missing some parsing info")
    return output 

#TAG::alignment score, miRanda::
def calculateAlignmentScore(query, reference):

    allLetters = ['A', 'T', 'C', 'G']
    queryGapPairs = ['-%s' % x for x in allLetters]
    referenceGapPairs = ['%s-' % x for x in allLetters]
    #TAG::genomic letter combinations,combinations,letters::
    guPairs = ['GT', 'TG']
    compPairs = ['AT', 'TA', 'CG', 'GC']
    misPairs = ['TT', 'TC', 'AG', 'AC', 'AA', 'CC', 'CT', 'CA', 'GG', 'GA'] 
    
    #calculate new score
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

    smallMMPos = []
    numM = numMM = numGU = 0
    gapShift = 0
    for i,pair in enumerate(matchPairs):
        if pair in queryGapPairs:
            gapShift -= 1
        elif pair in guPairs:
            numGU += 1
        elif pair in compPairs:
            numM += 1
        elif pair in misPairs:
            smallMMPos.append(i - gapShift)
            numMM += 1
        elif pair in referenceGapPairs:
            pass
            #want to keep this so I know if I missed any combinations
        elif 'N' in pair:
            return [0,0,0,0,0,0] #dont take N alignments
        else:
            print query
            print reference
            print pair
            raise NameError("COMBINATION NOT ACCOUNTED FOR!!!")

    score = (5*numM) + (2*numGU) - (3*numMM) - (8 * (numQGaps + numRGaps)) - (2 * (numExtensionsQ + numExtensionsR))
    return [score, numMM, numM, numGU, numQGaps + numRGaps, smallMMPos]

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

