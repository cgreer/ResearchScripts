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
parRun(30, 3, '/home/chrisgre/myLibs/cgAlignmentFlat.py', 'appendTInfoFlat', aFN, dFN)
parClean(aFN, 30)
print timer.split()

'''
print '...updating paired interactions: centered mismatches and center expression'
parRun(30, 5, '/home/chrisgre/scripts/endoClip/updateMismatchAndMiddleFlat.py', 'markCenterExpression', aFN, '/home/chrisgre/smallLibs/siRNA/degradome/wigsk1')
parClean(aFN, 30)
print timer.split()
parRun(30, 3, '/home/chrisgre/scripts/endoClip/updateMismatchAndMiddleFlat.py', 'markMismatchedPairs', aFN)
parClean(aFN, 30)
print timer.split()
'''

