import bioLibCG
import cgNexusFlat
from cgNexusFlat import Field
import cgDegPeak

class cgAlignment:

        sID = Field('int', None, 1)                                     
        tID = Field('int', None, 2)
        sStart = Field('int', None, 3)
        sEnd = Field('int', None, 4)
        tStart = Field('int', None, 5)
        tEnd = Field('int', None, 6)
        sLength = Field('int', None, 7)
        tLength = Field('int', None, 8)
        numMismatches = Field('int', None, 9)
        mismatchPositions = Field('intList', list(), 10)
        tTcc = Field('string', '.', 11)
        mismatchStatus = Field('boolList', [False, False, False], 12)
        centerExpression = Field('floatList', [0.0, 0.0, 0.0], 13)
        transcriptOverlap = Field('bool', False, 14)
        passed = Field('bool', False, 15)
        tELevel = Field('int', 0, 16)
        context = Field('string', '.', 17)
        repeat = Field('bool', False, 18)
        gScore = Field('int', 0, 19)
        targetSequence = Field('string', '.', 20)

def appendTInfo(aFN, degSmallFN, rn = None, tn = None):

        aNX = cgNexusFlat.Nexus(aFN, cgAlignment)
        aNX.load(['tID', 'tTcc'], [rn, tn])

        tID_aIDs = {}
        for aID in aNX.tID:
                tID_aIDs.setdefault(aNX.tID[aID], []).append(aID)

        
        f = open(degSmallFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = int(ls[0])
                for aID in tID_aIDs.get(id, list()):
                        aNX.tTcc[aID] = ls[1]

        aNX.save()                        

def appendTInfoFlat(aFN, dFN, rn = None, tn = None):

        aNX = cgNexusFlat.Nexus(aFN, cgAlignment)
        aNX.load(['tID', 'tTcc', 'transcriptOverlap', 'tELevel', 'context', 'repeat', 'targetSequence', 'gScore'], [rn, tn])

        dNX = cgNexusFlat.Nexus(dFN, cgDegPeak.Peak)
        dNX.load(['tOverlap', 'eLevel', 'tcc', 'context', 'repeatStatus', 'sequence', 'gScore'])

        tID_aIDs = {}
        for aID in aNX.tID:
                tID_aIDs.setdefault(aNX.tID[aID], []).append(aID)


        for dID in dNX.tcc: 
                
                for aID in tID_aIDs.get(dID, list()):
                        aNX.tTcc[aID] = dNX.tcc[dID]
                        aNX.tELevel[aID] = dNX.eLevel[dID]
                        aNX.transcriptOverlap[aID] = dNX.tOverlap[dID]
                        aNX.context[aID] = dNX.context[dID]
                        aNX.repeat[aID] = dNX.repeatStatus[dID]
                        aNX.gScore[aID] = dNX.gScore[dID]
                        aNX.targetSequence[aID] = dNX.sequence[dID]
        aNX.save()                        

def appendTranInfo(aFN, degSmallFN, rn = None, tn = None):

        aNX = cgNexusFlat.Nexus(aFN, cgAlignment)
        aNX.load(['tID', 'transcriptOverlap'], [rn, tn])

        # get id_overlap
        tID_tranVal = {}
        f = open(degSmallFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                val = int(ls[3])
                if val == 1: val = True
                if val == 0: val = False

                tID_tranVal[int(ls[0])] = val
        
        #update it.
        for aID in aNX.transcriptOverlap:
                tID = aNX.tID[aID]
                aNX.transcriptOverlap[aID] = tID_tranVal[tID]
                

        aNX.save()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])


