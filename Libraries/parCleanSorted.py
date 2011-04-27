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
        os.remove(fN)
        
        #remove the exit signals
        for eFN in bioLibCG.recurseDir(dirName, start = basename, end = 'exitSignal'):
                pass
                os.remove(eFN)

        #cat the SORTED range files.                
        rangeFiles = bioLibCG.recurseDir(dirName, start = basename, within = 'range')
        rangeFiles.sort(key = lambda x: int(x.split('.')[-2]))
        
        f = open(fN, 'w')
        for rFN in rangeFiles:
                subprocess.Popen(['cat', rFN], stdout = f).wait()
        f.close()

        #remove the range files
        for rFN in rangeFiles:
                os.remove(rFN)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
