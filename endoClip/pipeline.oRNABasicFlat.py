import bioLibCG
import siRnaPredictFlat as si
import sys
import updateContigsFlat 
import updateDuplicatesMultiTcc as udmt
import blankIDs

#this is the directory where the alignment file has been placed

'''Pipeline assumes that full alignment of small to target seqs has been done, but not truncated'''
'''Also assumes that the results files have been made...'''
oRNAFN = sys.argv[1]
peakFN = sys.argv[2]
seqFN = sys.argv[3]

timer = bioLibCG.cgTimer()
timer.start()

#initialize oRNA database
print 'make blank IDs'
blankIDs.blankIDs(seqFN, oRNAFN)
print 'sequence'
si.updateSequence(oRNAFN, seqFN)
print 'tcc'
si.updateTcc(oRNAFN, peakFN)
print 'entropy'
si.updateEntropy(oRNAFN)
si.updateELevel(oRNAFN, '/home/chrisgre/smallLibs/siRNA/small/wigsk50)

print 'contigs'
updateContigsFlat.updateTotalContig(oRNAFN)
updateContigsFlat.updateEndContig(oRNAFN)
print 'duplicates and mutlitcc'
udmt.updateSeqDuplicateMultiTcc(oRNAFN)
print timer.split()

