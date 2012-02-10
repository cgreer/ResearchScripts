import os
import bioLibCG
import subprocess

def parClean(fN):
        
        if '/' in fN:
                dirName = os.path.dirname(fN)
        else:
                dirName = os.environ['PWD']

        basename = os.path.basename(fN)

        #remove the original file
        print '..removing original file, if present'
        try:
            os.remove(fN)
        except OSError:
            pass
        
        #remove the exit signals
        print '..removing exit signals'
        for eFN in bioLibCG.recurseDir(dirName, start = basename, end = 'exitSignal'):
                os.remove(eFN)

        #cat the SORTED range files.                
        rangeFiles = bioLibCG.recurseDir(dirName, start = basename, within = 'range')
        rangeFiles.sort(key = lambda x: int(x.split('.')[-2]))
 

        print '..catting files together'
        f = open(fN, 'w')
        for rFN in rangeFiles:
                print '....', rFN
                subprocess.Popen(['cat', rFN], stdout = f).wait()
        f.close()

        #remove the range files
        print '..removing packets'
        for rFN in rangeFiles:
                print '....', rFN
                os.remove(rFN)

def parCleanSplit(fN):
        '''Remove the exit signals for a split continuation run'''
        if '/' in fN:
                dirName = os.path.dirname(fN)
        else:
                dirName = os.environ['PWD']

        basename = os.path.basename(fN)

        #remove the exit signals
        for eFN in bioLibCG.recurseDir(dirName, start = basename, end = 'exitSignal'):
                os.remove(eFN)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
