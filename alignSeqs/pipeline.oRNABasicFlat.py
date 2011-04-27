import bioLibCG
import siRnaPredictFlat as si
import sys
import updateContigsFlat 
import updateDuplicatesMultiTcc as udmt
import blankIDs

#this is the directory where the alignment file has been placed
runName = sys.argv[1]

'''Pipeline assumes that full alignment of small to target seqs has been done, but not truncated'''
'''Also assumes that the results files have been made...'''
oRNAFN = sys.argv[2]

peakFN = runName + '/oRNA.peaks'
seqFN = runName + '/oRNA.sequences'
tFN = runName + '/dRNA.results'

timer = bioLibCG.cgTimer()
timer.start()

#initialize oRNA database

print 'make blank IDs'
#blankIDs.blankIDs(seqFN, oRNAFN)
print 'sequence'
si.updateSequence(oRNAFN, seqFN)
print 'tcc'
si.updateTcc(oRNAFN, peakFN)
print 'entropy'
si.updateEntropy(oRNAFN)
#si.updateSmallExpression(oRNAFN, 'siPeaks.conf')
print 'contigs'
updateContigsFlat.updateTotalContig(oRNAFN)
updateContigsFlat.updateEndContig(oRNAFN)
print 'duplicates and mutlitcc'
udmt.updateSeqDuplicateMultiTcc(oRNAFN)
print timer.split()

