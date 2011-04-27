import bioLibCG
import cgAlign
import cgExit

def alignSeqs(seqsFN, dbName, wordSize, outFN, maxNumMismatches, sendExitSignal = False):
        maxNumMismatches = int(maxNumMismatches)
        sendExitSignal = bool(sendExitSignal)

        timer = bioLibCG.cgTimer()
        timer.start()
        #put seqs in cgSeq object, align
        wName = dbName + '.wDB'
        sName = dbName + '.sDB'
        wordSize = int(wordSize)
        
        #load dbs
        #print 'loading Sequence Database'
        sDB = cgAlign.loadSequenceDatabase(sName)
        print timer.split()
        #print 'loading Word Database'
        wDB = cgAlign.loadWordDatabase(wName)
        print timer.split()

        #align each seq
        f = open(seqsFN, 'r')
        fOut = open(outFN, 'w')
        for line in f:
                qSeq = cgAlign.cgSeq(line.strip().split('\t')[0], line.strip().split('\t')[1])
                
                #write out the alignments
                cgAlign.alignQuery(qSeq, wDB, sDB, wordSize, maxNumMismatches, fOut)
        
        f.close()
        fOut.close()
        
        print timer.split()
        if sendExitSignal:
                cgExit.sendExitSignal(seqsFN)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(alignSeqs, sys.argv)
