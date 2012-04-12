import bioLibCG
from cgNexus import Nexus
import cgDL
from cgAutoCast import autocast
from cgAutoKeyWord import autokey
from bioLibJA import subit
import os

def linkTargetIDs(oFN, oFF, aFN, aFF):

    oNX = Nexus(oFN, oFF)
    oNX.load(['filteredTargets'])

    #just give it some blanks
    if os.path.getsize(aFN) == 0:
        while oNX.nextID():
            oNX.filteredTargets = []
        oNX.save()
        return
    

    aNX = Nexus(aFN, aFF)
    aNX.load(['sID'])

    sID_aIDs = aNX.createMap('sID', 'id', False)
    
    for sID, aIDs in sID_aIDs.iteritems():
        oNX.id = sID
        oNX.filteredTargets = aIDs

    oNX.save()

def updateNumUFBS(oFN, oFF, aFN, aFF):

    oNX = Nexus(oFN, oFF)
    oNX.load(['filteredTargets', 'numUFBS'])

    #just give it some blanks
    if os.path.getsize(aFN) == 0:
        while oNX.nextID():
            oNX.numUFBS = 0
        oNX.save()
        return
    
    aNX = Nexus(aFN, aFF)
    aNX.load(['sigMask'])

    while oNX.nextID():
        sigMaskSet = set()
        for aID in oNX.filteredTargets:
            aNX.id = aID
            sigMaskSet.add(aNX.sigMask)
        oNX.numUFBS = len(sigMaskSet)

    oNX.save()

def updateAvgSS(dataFN, oFF, simDir, simBase, mm, numSims = 100):

    #get simulation information (# unique sims, num UFBS)
    fileNames = ['%s/simulation.%s/%s.%s' % (simDir, i, simBase, mm) for i in range(numSims)]
    sID_numSimUFBS = {}
    sID_simSeqs = {}
    for fN in fileNames:
        oNX = Nexus(fN, oFF)
        oNX.load(['numUFBS', 'sequence'])
        while oNX.nextID():
            if oNX.sequence in sID_simSeqs.get(oNX.id, set()):
                pass #dont count again
            else:
                sID_numSimUFBS.setdefault(oNX.id, []).append(oNX.numUFBS)
                sID_simSeqs.setdefault(oNX.id, set()).add(oNX.sequence)

    #update data based on sim info
    dataNX = Nexus(dataFN, oFF)
    dataNX.load(['avgNumSimSS', 'numUniqueSims', 'totalNumUFBSSim'])
    while dataNX.nextID():

        numUniqueSims = len(sID_numSimUFBS.get(dataNX.id, []))
        totalSimUFBS = sum(sID_numSimUFBS.get(dataNX.id, []))
        avgSimUFBS = totalSimUFBS/float(numUniqueSims) if numUniqueSims != 0 else -1
        dataNX.avgNumSimSS = avgSimUFBS
        dataNX.numUniqueSims = numUniqueSims
        dataNX.totalNumUFBSSim = totalSimUFBS 

    dataNX.save()

def updateSNR(oFN, oFF):

    NX = Nexus(oFN, oFF)
    NX.load(['avgNumSimSS', 'numUFBS', 'snr'])

    while NX.nextID():
        try:
            NX.snr = NX.numUFBS/NX.avgNumSimSS
        except ZeroDivisionError:
            NX.snr = NX.numUFBS/0.01

    NX.save()

def updateSimSeqsForUnique(oFN, oFF, seqFN):

    NX = Nexus(oFN, oFF)
    NX.load(['sequence'])

    id_seq = {}
    f = open(seqFN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        id_seq[int(ls[0])] = ls[1]
    f.close()
    

    while NX.nextID():
        NX.sequence = id_seq.get(NX.id, '.')

    NX.save()

def simFAtoSeqs(faFN, outFN):

    simSeqs = []
    simIDs = []
    f = open(faFN, 'r')
    for i,line in enumerate(f):
        if i %3 == 0:
            simIDs.append(line.strip()[1:])
        if i % 3 == 1:
            simSeqs.append(line.strip())
    f.close()

    fOut = open(outFN, 'w')
    for (id, seq) in zip(simIDs, simSeqs):
        fOut.write('%s\t%s\n' % (id, seq))

def createSubSequences(oID_sequence, frameLength):
 
    return dict( (oID, set(bioLibCG.returnFrames(oID_sequence[oID], frameLength))) for oID in oID_sequence)

def createSubMerOverlapMaps(oID_subMers):

    oID_connectedIDs = {}
    for oID, oSubset in oID_subMers.iteritems():
        oID_connectedIDs[oID] = set()
        for anotherID, anotherSet in oID_subMers.iteritems():
            if oID == anotherID: continue
            if oSubset.intersection(anotherSet):
                oID_connectedIDs[oID].add(anotherID)

    return oID_connectedIDs

def consolidateOverlapSets(oID_connectedIDs):
    [oID_connectedIDs[x].add(x) for x in oID_connectedIDs]
    allIDs = oID_connectedIDs.values()
    
    debugSet = set([2344, 1999, 2867, 2294])

    consolidatedSets = []
    for idSet in allIDs:
        addedToSet = False

        #consolidate if it shares constituents with another set
        for cSet in consolidatedSets:
            if idSet.intersection(cSet):
                [cSet.add(x) for x in cSet.union(idSet)]
                addedToSet = True
                break

        #if no consolidation took place, add it as its own set
        if not addedToSet:
            consolidatedSets.append(idSet)

    return consolidatedSets

def getSimilarORNASets(oID_sequence, frameLength):

    oID_subMers = createSubSequences(oID_sequence, frameLength)
    oID_connectedIDs = createSubMerOverlapMaps(oID_subMers)
    return consolidateOverlapSets(oID_connectedIDs)

@autocast
def testConsolidation(oFN, oFF, frameLength):

    dataNX = Nexus(oFN, oFF)
    dataNX.load(['sequence'])
    oID_sequence = dataNX.createMap('id', 'sequence')
    consolidatedSets = getSimilarORNASets(oID_sequence, frameLength)
    
    #check if all oIDs are in set
    allConsolidatedIDs = set()
    [allConsolidatedIDs.add(x) for theSet in consolidatedSets for x in theSet]    
    oIDsSet = set(oID_sequence.keys())
    print "DIFFERENCE"
    print oIDsSet.symmetric_difference(allConsolidatedIDs)

    #check Duplicates
     
    #print out sets to verify that they work 
    for oIDSet in consolidatedSets:
        print 
        print oIDSet
        for oID in oIDSet:
            print oID, oID_sequence[oID]

@autocast
def updateSimilarSiblings(oFN, oFF, frameLength):

    dataNX = Nexus(oFN, oFF)
    dataNX.load(['sequence', 'siblingSet'])
    oID_sequence = dataNX.createMap('id', 'sequence')
    consolidatedSets = getSimilarORNASets(oID_sequence, frameLength)

    for cSet in consolidatedSets:
        for oID in cSet:
            dataNX.id = oID
            dataNX.siblingSet = list(cSet)

    dataNX.save()

def testOverlaps(dataFN, oFF):

    dataNX = Nexus(dataFN, oFF)
    dataNX.load(['tcc'])

    #check for overlaps
    overlappingIDs = set()
    chrom_strand_range = {}
    while dataNX.nextID():
        chrom, strand, start, end = bioLibCG.tccSplit(dataNX.tcc)

        #check if overlap
        chrom_strand_range.setdefault(chrom, {}).setdefault(strand, set())
        overlap = False
        for i in range(start, end + 1):
            if i in chrom_strand_range[chrom][strand]:
                overlap = True
                break

        #tag or add these coordinates 
        if overlap:
            overlappingIDs.add(dataNX.id)
        else:
            for i in range(start, end + 1):
                chrom_strand_range[chrom][strand].add(i)

    print "THESE OVERLAP", overlappingIDs

def cleanForSNR(dataFN, oFF):
    dataNX = Nexus(dataFN, oFF)
    dataNX.load(['numUniqueSims', 'numUFBS', 'snrClean', 'siblingSet'])
    id_numUFBS = dataNX.createMap('id', 'numUFBS')
    id_siblingSet = dataNX.createMap('id', 'siblingSet')

    unusedSiblings = []
    for id, siblingSet in id_siblingSet.iteritems():
        if len(siblingSet) == 1: continue #NOTE: oRNA IDs are in their own sibling set
        numUFBS__id = [(id_numUFBS[x], x) for x in siblingSet]
        numUFBS__id.sort()
        numUFBS__id.pop() #take last one (one we're keeping) out of list 
        unusedIDs = [x[1] for x in numUFBS__id]
        unusedSiblings.extend(unusedIDs)

    #tag unclean oRNA
    while dataNX.nextID():
        if (dataNX.id in unusedSiblings) or (dataNX.numUniqueSims < 10):
            dataNX.snrClean = False
        else:
            dataNX.snrClean = True
    dataNX.save()

@autocast
def calculateTotalSNR(dataFN, oFF, mm, iSNRCutoff):
    mm = str(mm)
    
    dataNX = Nexus(dataFN, oFF)
    dataNX.load(['snr', 'numUFBS', 'numUniqueSims', 'totalNumUFBSSim']) 

    #Check SNR Cutoffs
    dataIDs = set()
    lowSNRORNA = set()
    while dataNX.nextID():
        dataIDs.add(dataNX.id)
        if dataNX.snr < iSNRCutoff:
            lowSNRORNA.add(dataNX.id)

    #NOTE sum of avgs != total avg
    #get total numUFBS for data and simulation total,num unique sims
    totalUFBSData = 0.0
    totalUFBSSim = 0.0
    totalPassingORNA = 0
    totalUniqueSims = 0
    while dataNX.nextID():
        if (dataNX.id in lowSNRORNA): continue
        totalUFBSData += dataNX.numUFBS
        totalPassingORNA += 1
        totalUFBSSim += dataNX.totalNumUFBSSim
        totalUniqueSims += dataNX.numUniqueSims

    totalAvgSimUFBS, totalAvgDataUFBS = 0.0, 0.0
    try:
        totalAvgSimUFBS = totalUFBSSim / float(totalUniqueSims)
        totalAvgDataUFBS = totalUFBSData / float(totalPassingORNA)
        totalSNR = totalAvgDataUFBS/totalAvgSimUFBS
        oS = [str(x) for x in [mm, iSNRCutoff, totalUFBSData, totalAvgSimUFBS, totalSNR, totalPassingORNA, float(totalUFBSData)/totalPassingORNA]]
    except ZeroDivisionError:
        oS = [str(x) for x in [mm, iSNRCutoff, totalUFBSData, totalAvgSimUFBS, "NA", totalPassingORNA, "NA"]]
        
    print '\t'.join(oS)

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

