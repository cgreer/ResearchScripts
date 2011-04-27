import bioLibCG
import siRnaPredict as si
import sys
import updateContigs
import updateDuplicatesMultiTcc as udmt

#this is the directory where the alignment file has been placed
runName = sys.argv[1]

'''Pipeline assumes that full alignment of small to target seqs has been done, but not truncated'''
'''Also assumes that the results files have been made...'''
oRNADir = runName + '/oRNA'

peakFN = 'oRNA.peaks'
seqFN = 'oRNA.sequences'
tFN = 'dRNA.results.e.t'

timer = bioLibCG.cgTimer()
timer.start()

#initialize oRNA database
print 'initializing oRNA'
si.updateID(oRNADir, peakFN)
print 'sequence'
si.updateSequence(oRNADir, seqFN)
print 'tcc'
si.updateTcc(oRNADir, peakFN)
print 'entropy'
si.updateEntropy(oRNADir)
#si.updateSmallExpression(oRNADir, 'siPeaks.conf')
print 'contigs'
updateContigs.updateTotalContig(oRNADir)
updateContigs.updateEndContig(oRNADir)
print 'duplicates and mutlitcc'
udmt.updateSeqDuplicateMultiTcc(oRNADir)
print timer.split()

