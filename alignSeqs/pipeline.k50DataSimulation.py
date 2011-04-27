import bioLibCG
import sys
import blankIDs
from parRun import parRun
from checkExit import parClean

#this is the directory where the alignment file has been placed
runName = sys.argv[1]

'''Pipeline assumes that full alignment of small to target seqs has been done, but not truncated'''
'''Also assumes that the results files have been made...'''
oFN = runName + '/oRNA.data'
aFN = runName + '/all.aligned.0.5.tabs'

seqFN = runName + '/oRNA.sequences'
tFN = runName + '/dRNA.results'

timer = bioLibCG.cgTimer()
timer.start()

#initialize oRNA
blankIDs.blankIDs(seqFN, oFN)

#initializing alignments
print '...appending T Info'
parRun(40, 3, '/home/chrisgre/scripts/alignSeqs/cgAlignmentFlat.py', 'appendTInfo', aFN, tFN)
parClean(aFN, 40)
print '...appending Tran Info'
parRun(40, 3, '/home/chrisgre/scripts/alignSeqs/cgAlignmentFlat.py', 'appendTranInfo', aFN, tFN)
parClean(aFN, 40)
print timer.split()


print '...updating paired interactions: centered mismatches and center expression'
parRun(30, 4, '/home/chrisgre/scripts/alignSeqs/updateMismatchAndMiddleFlat.py', 'markCenterExpression', aFN, '/home/chrisgre/smallLibs/siRNA/degradome/wigsk50')
parClean(aFN, 30)
print timer.split()
parRun(50, 3, '/home/chrisgre/scripts/alignSeqs/updateMismatchAndMiddleFlat.py', 'markMismatchedPairs', aFN)
parClean(aFN, 50)
print timer.split()

print 'Linking oRNA to Targets'
#update the initial targets for small rnas based off of only local alignment
#parRun(30, 4, '/home/chrisgre/scripts/alignSeqs/siRnaPredictFlat.py', 'updateTargetIDs', oFN, aFN)
#parClean(oFN, 30)
print timer.split(), 'done', runName


