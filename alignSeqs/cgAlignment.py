import bioLibCG
import cgDB

class cgAlignment:

        id = cgDB.Field('int', None)
        sID = cgDB.Field('int', None)
        tID = cgDB.Field('int', None)
        sStart = cgDB.Field('int', None)
        sEnd = cgDB.Field('int', None)
        tStart = cgDB.Field('int', None)
        tEnd = cgDB.Field('int', None)
        sLength = cgDB.Field('int', None)
        tLength = cgDB.Field('int', None)
        numMismatches = cgDB.Field('int', None)
        mismatchPositions = cgDB.Field('intList', [])
        tTcc = cgDB.Field('string', None)
        mismatchStatus = cgDB.Field('boolList', [False, False, False])
        centerExpression = cgDB.Field('floatList', [0.0, 0.0, 0.0])
        transcriptOverlap = cgDB.Field('bool', False)

        def __init__(self, id):
                self.id = id
                        
def pPrint(aDir):

        aDC = cgDB.dataController(aDir, cgAlignment)
        id_alignment = aDC.load()
        
        for alignment in id_alignment.values():
                attName_att = alignment.__dict__
                attVals = [attName_att['id'], attName_att['tID'], attName_att['tTcc']]
                attVals = [str(x) for x in attVals]
                print '\t'.join(attVals)

def prettyPrint(alignment):

        print alignment.id, alignment.sID, alignment.tID, alignment.centerExpression, alignment.mismatchStatus, alignment.numMismatches, alignment.transcriptOverlap

def loadAlignments(aDir, alignmentFN):
        
        aDC = cgDB.dataController(aDir, cgAlignment)

        id_alignment = {}
        i = 0
        f = open(alignmentFN, 'r')
        for line in f:
                ls = line.strip().split(' ')
                
                a = cgAlignment(i)
                
                a.sID, a.tID = int(ls[0]), int(ls[1])
                a.sStart, a.sEnd = int(ls[2]), int(ls[3])
                a.tStart, a.tEnd = int(ls[4]), int(ls[5])
                a.sLength, a.tLength = int(ls[6]), int(ls[7])
                a.numMismatches = int(ls[8])
                
                try:
                        a.mismatchPositions = [int(x) for x in ls[9].split(',')]
                except IndexError:
                        a.mismatchPositions = []

                id_alignment[i] = a

                i += 1
        
        aDC.commit(id_alignment)

def loadAlignments2(aDir, alignmentFN):
        '''Just added IDs at the beginning to parallel things'''
        
        aDC = cgDB.dataController(aDir, cgAlignment)

        id_alignment = {}
        f = open(alignmentFN, 'r')
        for line in f:
                ls = line.strip().split(' ')
                id = int(ls[0]) 
                a = cgAlignment(id)
                
                a.sID, a.tID = int(ls[1]), int(ls[2])
                a.sStart, a.sEnd = int(ls[3]), int(ls[4])
                a.tStart, a.tEnd = int(ls[5]), int(ls[6])
                a.sLength, a.tLength = int(ls[7]), int(ls[8])
                a.numMismatches = int(ls[9])
                
                try:
                        a.mismatchPositions = [int(x) for x in ls[10].split(',')]
                except IndexError:
                        a.mismatchPositions = []

                id_alignment[id] = a

        aDC.commit(id_alignment)

def appendTInfo(aDir, degSmallFN):

        aDC = cgDB.dataController(aDir, cgAlignment)
        id_alignment = aDC.load()

        tID_alignments = {}
        for alignment in id_alignment.values():
                tID_alignments.setdefault(alignment.tID, []).append(alignment)

        
        f = open(degSmallFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = int(ls[0])
                for a in tID_alignments.get(id, []):
                        a.tTcc = ls[1]
                        if id == 192002:
                                print a.tTcc, a.id

        aDC.commit(id_alignment)

def appendTranInfo(aDir, degSmallFN):

        aDC = cgDB.dataController(aDir, cgAlignment)
        id_alignment = aDC.load()

        tID_tranVal = {}
        f = open(degSmallFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                val = int(ls[3])
                if val == 1: val = True
                if val == 0: val = False

                tID_tranVal[int(ls[0])] = val
        
        for alignment in id_alignment.values():
                transcriptOverlap = tID_tranVal[alignment.tID]
                alignment.transcriptOverlap = transcriptOverlap
                

        aDC.commit(id_alignment)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(loadAlignments, sys.argv)


