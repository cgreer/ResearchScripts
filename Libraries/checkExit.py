import bioLibCG
import os
import time
import subprocess
import parCleanSorted

def checkExit(fN, numPackets):
        numPackets = int(numPackets)
        dirName = os.path.dirname(fN)
        baseName = os.path.basename(fN)
        if os.environ['PWD'] not in dirName:
                dirName = os.path.dirname(os.environ['PWD'] + '/' + fN)
 
        sleepTime = 20
        iteration = 1
        while True:

                time.sleep(sleepTime)
                exitSignals = bioLibCG.recurseDir(dirName, start = baseName , end = 'exitSignal')
                print 'waiting...', str(len(exitSignals)), '/', str(numPackets), '%s' % bioLibCG.prettyTime(sleepTime * iteration), fN 
                if len(exitSignals) == numPackets:
                        print 'Jobs all finished!'
                        break
                iteration += 1

def parClean(fN, numPackets):
        
        #wait until done
        checkExit(fN, numPackets)

        #run cleanup script
        home = os.environ['HOME']
        parCleanSorted.parClean(fN)


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
