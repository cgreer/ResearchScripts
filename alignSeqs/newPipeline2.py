import bioLibCG
import siRnaPredict as si
import cgAlignment
import updateMismatchAndMiddle as mm
import sys

#this is the directory where the alignment file has been placed
runName = sys.argv[1]
'''Pipeline assumes that full alignment of small to target seqs has been done, but not truncated'''
'''Also assumes that the results files have been made...'''
oRNADir = runName + '/oRNA'
aDir = runName + '/aDir'

minMismatches = 0
maxMismatches = 5

alignFN = runName + '/all.aligned'
alignTrunFN = runName + '/all.aligned' + '.' + str(minMismatches) + '.' + str(maxMismatches)
peakFN = 'oRNA.results'
tFN = 'dRNA.results.e.t'

timer = bioLibCG.cgTimer()
timer.start()

#initialize oRNA database
print 'initializing oRNA'
si.updateID(oRNADir, peakFN)
print timer.split()

print 'truncating alignments'
#Truncate alignment file for specified mismatches
#si.truncateAlignments(alignFN, minMismatches, maxMismatches, alignTrunFN)

print 'loading alignments'
#initialize alignments
#cgAlignment.loadAlignments2(aDir, alignTrunFN)
#cgAlignment.appendTInfo(aDir, tFN)
#cgAlignment.appendTranInfo(aDir, tFN)
print timer.split()

print 'updating initial targets'
#update the initial targets for small rnas based off of only local alignment
#si.updateTargetIDs(oRNADir, aDir)
print timer.split()

print 'updating paired interactions: centered mismatches and center expression'
mm.markCenterExpression(aDir, 'siDegradome.conf')
mm.markMismatchedPairs(aDir)
print timer.split()
