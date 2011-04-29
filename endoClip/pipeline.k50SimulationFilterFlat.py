import bioLibCG
import sys
from parRun import parRun
from checkExit import parClean
import subprocess
import copyIDs

#fN and Init
oFN = sys.argv[1]
aFN = sys.argv[2] 
seqs = sys.argv[3]
aFilteredFN = aFN  #+ '.filtered'

timer = bioLibCG.cgTimer()
timer.start()

#filter targets
print 'filtering targets'
#parRun(30, 3, '/home/chrisgre/scripts/endoClip/filteringFlat.py', 'filterTargetsInPlace', aFN, 'True', '1', '1', '.60')
#parClean(aFN, 30)
print timer.split()

#make db smaller
print 'truncating db'
#subprocess.Popen(['/home/chrisgre/scripts/endoClip/truncate.filterAlignments.sh', aFN, aFilteredFN]).wait()
print timer.split()

#make copy of ids
copyIDs.copyIDs(seqs, oFN)

#link oRNA to filteredTargets
print 'linking oRNA to filtered targets'
parRun(1, 'LOCAL', '/home/chrisgre/scripts/endoClip/siRnaPredictFlat.py', 'updateTargetIDsFiltered', oFN, aFilteredFN)
#parClean(oFN, 30)
print timer.split()


