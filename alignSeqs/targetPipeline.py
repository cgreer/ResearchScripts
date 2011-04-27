import bioLibCG
import siRnaPredict as si
import updateMismatchAndMiddle as mm
import sys

runName = 'test/test'
runName = sys.argv[1]
'''Pipeline assumes that full alignment of small to target seqs has been done, but not truncated'''
'''Also assumes that the results files have been made...'''

minMismatches = 0
maxMismatches = 5

alignFN = runName + '.all.aligned.nodups'
alignTrunFN = runName + '.all.aligned.nodups' + '.' + str(minMismatches) + '.' + str(maxMismatches)
sFN = runName + '.small.results'
tFN = runName + '.degradome.results'


#Get expression of the small RNAs and the targets
print 'updating expression levels'
#si.updateSmallExpression(sFN, 'siPeaks.conf', 1, 3, sFN) 
#si.updateSmallExpression(tFN, 'siDegradome.conf', 1, 2, tFN)

print 'truncating alignments'
#Truncate alignment file for specified mismatches
si.truncateAlignments(alignFN, minMismatches, maxMismatches, alignTrunFN)


print 'updating initial targets'
#update the initial targets for small rnas based off of only local alignment
si.updateTargetIDs(sFN, alignTrunFN, 0, 4, sFN)

print 'update microrna/transcript overlaps'
#update microRNA biogenesis for small and transcript
#si.updateMicroRNAOverlap(sFN, 'mirBaseHumanNew.gff.Tcc', 1, 5, sFN)
#si.transcriptSetOverlap(sFN, False, 1, 6, sFN) 
#si.transcriptSetOverlap(tFN, True, 1, 3, tFN)


print 'updating paired interactions: centered mismatches and center expression'
mm.markCenterExpression(sFN, tFN, alignTrunFN, 'siDegradome.conf', runName + '.pair.center.data')
mm.markMismatchedPairs(sFN, alignTrunFN, runName + '.pair.mismatch.data')

