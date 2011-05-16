import bioLibCG
import sys
from parRun import parRun
from checkExit import parClean
import subprocess

#fN and Init
oFN = sys.argv[1]
aFN = sys.argv[2] #this is the filtered ONE! 

timer = bioLibCG.cgTimer()
timer.start()


#link oRNA to filteredTargets
print 'linking oRNA to filtered targets'
parRun(1, 'LOCAL', '/home/chrisgre/scripts/endoClip/siRnaPredictFlat.py', 'updateTargetIDsFiltered', oFN, aFN)
print timer.split()


