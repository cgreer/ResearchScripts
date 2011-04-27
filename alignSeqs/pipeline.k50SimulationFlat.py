import bioLibCG
import sys
import blankIDs
from parRun import parRun
from checkExit import parClean


#files and init
oFN = sys.argv[1]
aFN = sys.argv[2] 

seqFN = sys.argv[3] 
tFN = '/home/chrisgre/scripts/alignSeqs/dRNA.results.updated'

timer = bioLibCG.cgTimer()
timer.start()


#initialize oRNA
blankIDs.blankIDs(seqFN, oFN)

#initializing alignments
print '...appending T Info'
parRun(30, 3, '/home/chrisgre/scripts/alignSeqs/cgAlignmentFlat.py', 'appendTInfo', aFN, tFN)
parClean(aFN, 30)
print '...appending Tran Info'
parRun(30, 3, '/home/chrisgre/scripts/alignSeqs/cgAlignmentFlat.py', 'appendTranInfo', aFN, tFN)
parClean(aFN, 30)
print timer.split()


print '...updating paired interactions: centered mismatches and center expression'
parRun(30, 5, '/home/chrisgre/scripts/alignSeqs/updateMismatchAndMiddleFlat.py', 'markCenterExpression', aFN, '/home/chrisgre/smallLibs/siRNA/degradome/wigsk50')
parClean(aFN, 30)
print timer.split()
parRun(30, 3, '/home/chrisgre/scripts/alignSeqs/updateMismatchAndMiddleFlat.py', 'markMismatchedPairs', aFN)
parClean(aFN, 30)
print timer.split()

print 'Linking oRNA to Targets'
'''will do this when I filter'''
#update the initial targets for small rnas based off of only local alignment
#parRun(30, 3, '/home/chrisgre/scripts/alignSeqs/siRnaPredictFlat.py', 'updateTargetIDs', oFN, aFN)
#parClean(oFN, 30)
print timer.split(), 'done', aFN


