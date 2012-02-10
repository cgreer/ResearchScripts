import bioLibCG
import sys

#!!!NOTE, the GT order is important.  The target db must contain the 3'->5'(comp) sequences.
#the queries are 5'->3'(normal)
rnaPairsGU = {'A': ['A', 'G'],
              'T': ['T'],
              'C': ['C', 'T'],
              'G': ['G'],
              'N':[]
             }

#for plain matches
#rnaPairsGU = {'A': ['A'],
              #'T': ['T'],
              #'C': ['C'],
              #'G': ['G'],
              #'N':[]
             #}

class cgAlignment:
        def __init__(self):
                self.target = None
                self.query = None
                self.qStart = None
                self.qEnd = None
                self.tStart = None
                self.tEnd = None
                self.numMismatches = None
                self.qID = None
                self.tID = None
                self.qLength = None
                self.tLength = None
                self.mismatchPositions = None

        def alignmentOutput(self):
                misMP = ','.join(self.mismatchPositions)
                infoString = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (self.qID, self.tID, self.qStart, self.qEnd, self.tStart, self.tEnd, self.qLength, self.tLength, self.numMismatches, misMP)
                return infoString 


class cgSeq:

        def __init__(self, id, sequence):
                self.id = id
                self.sequence = sequence

def createSequenceDatabase(seqList):
        seqDict = {} #: id: sequence
        for seq in seqList:
                seqDict[seq.id] = seq.sequence

        return seqDict

def writeSequenceDatabase(seqDict, runName):
        f = open(runName + '.sDB', 'w')
        for seqID in seqDict:
                f.write('%s\t%s\n' % (seqID, seqDict[seqID]))
        f.close()

def loadSequenceDatabase(fN):
        f = open(fN, 'r')

        seqDict = {} # id: sequence
        for line in f:
                ls = line.strip().split('\t')
                seqDict[int(ls[0])] = ls[1]
        f.close()

        return seqDict

def createWordDatabase(seqList, wordSize):
        #targetList should be a list of cgSeqs

        wordDict = {} # word: [id, word start, word end]
        for tSeq in seqList:
                for i,word in enumerate(bioLibCG.returnFrames(tSeq.sequence, wordSize)):
                        try:
                                wordDict[word].append([tSeq.id, i, i + wordSize - 1])
                        except KeyError:
                                wordDict[word] = [[tSeq.id, i, i + wordSize - 1]]
        return wordDict

def writeWordDatabase(wordDict, runName):
        f = open(runName + '.wDB', 'w')
        for word in wordDict:
                f.write(word + '\n')
                for wordPosition in wordDict[word]:
                        f.write('%s\t%s\t%s\n' % (wordPosition[0], wordPosition[1], wordPosition[2]))
        f.close()

def loadWordDatabase(fN):

        f = open(fN, 'r')
        
        wordDict = {}
        currentWord = None
        currentPositions = []
        for line in f:
                if not '\t' in line: #word line, start new word, end old one
                        if not currentWord == None:
                                wordDict[currentWord] = currentPositions
                                currentWord = line.strip()
                                currentPositions = []
                        else:
                                currentWord = line.strip()
                                currentPositions = []
                else: #word position, add to current word
                        currentPositions.append([int(line.strip().split('\t')[0]), int(line.strip().split('\t')[1]), int(line.strip().split('\t')[2])])
        f.close()

        #add final word to dict
        wordDict[currentWord] = currentPositions

        return wordDict

def alignSeqs(qSeq, tSeq, qWordPos, tWordPos, maxNumMismatches):
        #qWordPos is list [a,b] where a is start of word, b is end
        
        pairs = rnaPairsGU        
        query = qSeq.sequence
        target = tSeq.sequence
        alignment = cgAlignment()
        alignment.qID = qSeq.id
        alignment.tID = tSeq.id
        alignment.numMismatches = 0
        alignment.query = query
        alignment.target = target
        alignment.qLength = len(query)
        alignment.tLength = len(target)
        alignment.mismatchPositions = []

        #extend to the left
        qLeftRange = range(0, qWordPos[0])
        tLeftRange = range(0, tWordPos[0])
        
        i = qWordPos[0]
        j = tWordPos[0]
        prevI = i
        prevJ = j

        while True:
                try:
                        i = qLeftRange.pop()
                        j = tLeftRange.pop()
                        
                        #print query[i], target[j]
                        if query[i] not in pairs[target[j]]:
                                alignment.numMismatches += 1
                                if alignment.numMismatches > maxNumMismatches:
                                        return None
                                alignment.mismatchPositions.append(str(i))
                        
                        #This takes care of off by one popping exception...
                        prevI = i
                        prevJ = j

                except IndexError:
                        alignment.qStart = prevI
                        alignment.tStart = prevJ
                        break

        #make mismatchPosition in order
        alignment.mismatchPositions.reverse()
        
        #extend to the right
        qRightRange = range(qWordPos[1] + 1, len(query))
        qRightRange.reverse()
        tRightRange = range(tWordPos[1] + 1, len(target))
        tRightRange.reverse()


        i = qWordPos[1]
        j = tWordPos[1]
        prevI = i
        prevJ = j

        while True:

                try:
                        i = qRightRange.pop()
                        j = tRightRange.pop()
                        
                        if query[i] not in pairs[target[j]]:
                                alignment.numMismatches += 1
                                if alignment.numMismatches > maxNumMismatches:
                                        return None
                                alignment.mismatchPositions.append(str(i))
                        
                        prevI = i
                        prevJ = j

                except IndexError:
                        alignment.qEnd = prevI
                        alignment.tEnd = prevJ
                        break
        
        #Calculate residual mismatches from overhangs
        residual = len(query) - (1 + alignment.qEnd - alignment.qStart)
        alignment.numMismatches += residual
        
        if alignment.numMismatches > maxNumMismatches:
                return None
        else:
                return alignment

def alignQuery(qSeq, wordDatabase, sequenceDatabase, wordSize, maxNumMismatches, fOut):
        #wordDataBase : [id, wordStart, wordEnd]F
        #find all words in sequence --> align if word in database

        for i,word in enumerate(bioLibCG.returnFrames(qSeq.sequence, wordSize)):
                if word in wordDatabase:
                        for aInfo in wordDatabase[word]:
                                qWordPos = [i, i + wordSize - 1]
                                tSeq = cgSeq(aInfo[0], sequenceDatabase[aInfo[0]])
                                tWordPos = [aInfo[1], aInfo[2]]
                                newAlignment = alignSeqs(qSeq, tSeq, qWordPos, tWordPos, maxNumMismatches)
                                if newAlignment:
                                        fOut.write(newAlignment.alignmentOutput() + '\n')

                
        
if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(alignSeqs, sys.argv)





