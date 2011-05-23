import subprocess
import bioLibCG
import os


def splitRun(baseFN, memoryAmount, scriptName, *args):

        if '/' in baseFN:
                dirName = os.path.dirname(baseFN)
        else:
                dirName = os.environ['PWD']

        basename = os.path.basename(baseFN)

        #cat the SORTED range files.                
        rangeFiles = [x for x in bioLibCG.recurseDir(dirName, start = basename, within = 'range') if 'exitSignal' not in x]
        rangeFiles.sort(key = lambda x: int(x.split('.')[-2]))

        for fN in rangeFiles:

                #specific the correct qJob with correct memory
                qJobX = '%s/exec/qJobX%s.sh' % (os.environ['HOME'], memoryAmount)
                qDo = '%s/exec/qDo.sh' % (os.environ['HOME'])

                #construct command to pass
                if memoryAmount == 'LOCAL':
                        com = [qDo, scriptName]
                else:
                        com = [qJobX, qDo, scriptName]
                
                #append script arguments including SPLIT FN
                for arg in args:
                        if arg == 'splitFN':
                                com.append(fN)
                        else:                                
                                com.append(arg)
                        
                #run each job
                print com
                #subprocess.Popen(com).wait()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(splitRun, sys.argv)                
