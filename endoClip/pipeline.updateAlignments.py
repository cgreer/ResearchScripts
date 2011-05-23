import bioLibCG
import sys
import blankIDs
from parRun import parRun
from checkExit import parClean


#files and init
aFN = sys.argv[1] 
dFN = sys.argv[2]

timer = bioLibCG.cgTimer()
timer.start()


#initializing alignments
print '...appending T Info'
parRun(50, 3, '/home/chrisgre/myLibs/cgAlignmentFlat.py', 'appendTInfoFlat', aFN, dFN)
parClean(aFN, 50)
print timer.split()

print '...updating paired interactions: centered mismatches and center expression'
parRun(40, 5, '/home/chrisgre/scripts/endoClip/updateMismatchAndMiddleFlat.py', 'markCenterExpression', aFN, '/home/chrisgre/smallLibs/siRNA/degradome/wigsk50')
parClean(aFN, 40)
print timer.split()

parRun(50, 3, '/home/chrisgre/scripts/endoClip/updateMismatchAndMiddleFlat.py', 'markMismatchedPairs', aFN)
parClean(aFN, 50)
print timer.split()
