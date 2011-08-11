import bioLibCG
import os
import time
import subprocess
import parCleanSorted
from progressbar import ProgressBar,SimpleProgress,Percentage,Bar,Timer


def checkExit(fN, numPackets):
        numPackets = int(numPackets)
        dirName = os.path.dirname(fN)
        baseName = os.path.basename(fN)
        if os.environ['PWD'] not in dirName:
                dirName = os.path.dirname(os.environ['PWD'] + '/' + fN)
 
        sleepTime = 1
        pbar = ProgressBar(widgets=['  ', SimpleProgress(), ' ', Timer(), ' ', Bar()], maxval=numPackets).start()
        while True:
                time.sleep(sleepTime)
                exitSignals = bioLibCG.recurseDir(dirName, start = baseName , end = 'exitSignal')
                pbar.update(len(exitSignals))                
                if len(exitSignals) == numPackets:
                        pbar.finish()
                        print 'Jobs all finished!'
                        break

def parClean(fN, numPackets):
        
        #wait until done
        checkExit(fN, numPackets)

        #run cleanup script
        print '  cleaning up and merging files'
        parCleanSorted.parClean(fN)

def parCleanSplit(fN, numPackets):

        #wait until done
        checkExit(fN, numPackets)

        #run cleanup script
        print '  cleaning up for split run'
        parCleanSorted.parCleanSplit(fN)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
