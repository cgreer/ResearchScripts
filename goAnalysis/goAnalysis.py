import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast
from geneProperty import geneProperty
import cgDictHelp
import random
import numpy as np
from scipy.special import ndtr

def getPValInfoNormal(observation, samplingDist):

    sDist = np.array(samplingDist)
    theStd, theMean = sDist.std(), sDist.mean()
    zScore = (float(observation) - theMean)/theStd
    thePVal = ndtr(-zScore)
    return (theStd, theMean, zScore, thePVal)

def getPValInfo(observation, samplingDist):

    allNumbers = samplingDist
    allNumbers.append(observation)

    allNumbers.sort()
    numGreater = len([x for x in allNumbers if x >= observation])
    ePV = float(numGreater)/len(allNumbers)

    sDist = np.array(samplingDist)
    theStd, theMean = sDist.std(), sDist.mean()
    zScore = (float(observation) - theMean)/theStd
    return (theStd, theMean, zScore, ePV)

def getPVals(goInput, simCountFN, geneGoTerms): 

    #load data 
    gName_goTerms = cgDictHelp.dictFromColumns(geneGoTerms, (0, 1), ('string', 'stringList'))
    gName_backgroundSet = cgDictHelp.dictFromColumns(goInput, (0, 1), ('string', 'stringList'))
    goTerm_simCounts = cgDictHelp.dictFromColumns(simCountFN, (0, 1), ('string', 'intList'))

    #get data counts (observations)
    goTerm_count = dict( (term, 0) for term in goTerm_simCounts )
    updateGoCounts(goTerm_count, gName_goTerms, gName_backgroundSet.keys()) 

    #get pVals
    for goTerm, simCounts in goTerm_simCounts.iteritems(): 
        actualObservation = goTerm_count[goTerm]
        print '%s\t%s\t%s' % (goTerm, actualObservation, '\t'.join([str(x) for x in getPValInfo(actualObservation, simCounts)]))
        
def updateGoCounts(goTerm_count, gName_goTerms, geneList):

    for gName in geneList:
        releventTerms = [term for term in gName_goTerms[gName] if term in goTerm_count]
        for term in releventTerms:
            goTerm_count[term] += 1

@autocast
def simulateCounts(goInput, geneGoTerms, outFN, numSims = 1000):

    #collect all terms
    gName_goTerms = cgDictHelp.dictFromColumns(geneGoTerms, (0, 1), ('string', 'stringList'))

    #collect data set
    gName_backgroundSet = cgDictHelp.dictFromColumns(goInput, (0, 1), ('string', 'stringList'))
    
    #collect unique go terms in data set
    dataGoTerms = []
    for gName in gName_backgroundSet:
        dataGoTerms.extend(gName_goTerms[gName])
    dataGoTerms = list(set(dataGoTerms))

    #this will contain count for each sim performed
    goTerm_simCounts = dict( (term, []) for term in dataGoTerms )

    for i in range(1,numSims):
        #init counts
        goTerm_count = dict( (term, 0) for term in dataGoTerms)
        pairedControl = []
        for gName in gName_backgroundSet:
            pairedControl.append(random.choice(gName_backgroundSet[gName]))

        updateGoCounts(goTerm_count, gName_goTerms, pairedControl)
        for goTerm, count in goTerm_count.items():
            goTerm_simCounts[goTerm].append(count)

    f = open(outFN, 'w')
    for term in sorted(goTerm_simCounts.keys()):
        f.write('%s\t%s\n' % (term, ','.join([str(x) for x in goTerm_simCounts[term]])) )
    f.close()
    
def pickBackgroundGenes(goSet, geneSet, goInput, geneGoTerms):

    NX = cgNexusFlat.Nexus(geneSet, geneProperty)
    NX.load(['geneName', 'geneLength', 'geneGCContent'])
    
    #inverse on gName
    geneName_id = NX.iDictionary('geneName')
    
    #gene go terms to filter genes that dont have go terms
    gName_goTerms = cgDictHelp.dictFromColumns(geneGoTerms, (0, 1), ('string', 'stringList'))

    #first get length, GC ranges for go set
    usedGenes = []
    gName__lRange__GCRange = []
    f = open(goSet, 'r')
    for line in f:
        ls = line.strip().split('\t')
        gName = ls[0]
        usedGenes.append(gName)

        if gName not in geneName_id:
            print gName, 'not in dataset!'
            continue
        
        if gName not in gName_goTerms:
            print gName, 'has no entry in go terms!'
            continue

        nID = geneName_id[gName]
        gLength = NX.geneLength[nID]
        gGC = NX.geneGCContent[nID]
        
        gName__lRange__GCRange.append((gName, (gLength - .05*gLength, gLength + .05*gLength), (gGC - .05*gGC, gGC + .05*gGC)))

    f.close()

    #go through background set and pick acceptable genes 
    gName_bSet = {}
    for gName, lRange, GCRange in gName__lRange__GCRange:
        for id in NX.ids:
            bGeneName = NX.geneName[id]
            if bGeneName in usedGenes: continue
            if bGeneName not in gName_goTerms: continue
            theLength, theGC = NX.geneLength[id], NX.geneGCContent[id]

            if (GCRange[0] <= theGC <= GCRange[1]) and (lRange[0] <= theLength <= lRange[1]):
                gName_bSet.setdefault(gName, []).append(bGeneName)

    #write it to background file
    fOut = open(goInput, 'w')
    for gName, bSet in gName_bSet.iteritems():
        print gName, '# background genes', len(bSet)
        outString = [gName]
        outString.append(','.join(bSet))
        outString = '\t'.join([str(x) for x in outString])
        fOut.write(outString + '\n')
    fOut.close()
            
if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

