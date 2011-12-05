import bioLibCG
import cgAlign
from cgAutoCast import autocast

def createDatabases(targetsFN, wordSize, runName, hasIDs = False):
        wordSize = int(wordSize)
        hasIDs = (hasIDs == "True")
        print 'using IDs', hasIDs
        #make sequence list out of targets, make db, write to file
        f = open(targetsFN, 'r')
        seqList = []
        print 'obtaining sequences'
        i = 0
        for line in f:
                if hasIDs:
                    theID, seq = line.strip().split('\t')
                else:
                    theID, seq = i, line.strip()
                seqList.append(cgAlign.cgSeq(theID, seq))
                i += 1
        f.close()

        print 'making word db'
        wordDB = cgAlign.createWordDatabase(seqList, wordSize)
        cgAlign.writeWordDatabase(wordDB, runName)

        print 'making seq db'
        seqDB = cgAlign.createSequenceDatabase(seqList)
        cgAlign.writeSequenceDatabase(seqDB, runName)

        
if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(createDatabases, sys.argv)

