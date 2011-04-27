import bioLibCG
import cgAlign

def createDatabases(targetsFN, wordSize, runName):
        wordSize = int(wordSize)

        #make sequence list out of targets, make db, write to file
        f = open(targetsFN, 'r')
        seqList = []
        print 'obtaining sequences'
        i = 0
        for line in f:
                seqList.append(cgAlign.cgSeq(i, line.strip()))
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

