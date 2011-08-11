import bioLibCG
import sys
import blankIDs
from parRun import parRun
from checkExit import parClean
from checkExit import parCleanSplit
from splitRun import splitRun

#files and init
aFN = sys.argv[1] 
dFN = sys.argv[2]

timer = bioLibCG.cgTimer()
timer.start()


#initializing alignments
print '...appending T Info'
parRun(45, 3, '/home/chrisgre/myLibs/cgAlignmentFlat.py', 'appendTInfoFlat', aFN, dFN)
parCleanSplit(aFN, 45)
print timer.split()

print '...updating paired interactions:  center expression'
splitRun(aFN, 5, '/home/chrisgre/scripts/endoClip/updateMismatchAndMiddleFlat.py', 'markCenterExpression', 'splitFN', '/home/chrisgre/smallLibs/siRNA/degradome/wigsk50')
parCleanSplit(aFN, 45)
print timer.split()

print '...updating paired interactions:  center mismatch'
splitRun(aFN, 3, '/home/chrisgre/scripts/endoClip/updateMismatchAndMiddleFlat.py', 'markMismatchedPairs', 'splitFN')
parClean(aFN, 45)
print timer.split()
