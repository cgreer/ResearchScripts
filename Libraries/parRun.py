import subprocess
import bioLibCG
import os

def parRun(numParts, memoryAmount, scriptName, *args):
        numParts = int(numParts)
        noOut = open('/dev/null', 'w')
       
        print '  submitting Jobs'
        for i in xrange(1, numParts + 1):
               
                #specific the correct qJob with correct memory
                qJobX = '%s/exec/qJobX%s.sh' % (os.environ['HOME'], memoryAmount)
                qDo = '%s/exec/qDo.sh' % (os.environ['HOME'])

                #construct command to pass
                if memoryAmount == 'LOCAL':
                        com = [qDo, scriptName]
                else:
                        com = [qJobX, qDo, scriptName]
                
                #append script arguments
                for arg in args:
                        com.append(arg)
                        
                #append packet running info
                if numParts != 1:
                        com.append(str(i))
                        com.append(str(numParts))
                
                #run each job
                subprocess.Popen(com, stdout=noOut).wait()
        noOut.close()

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])



