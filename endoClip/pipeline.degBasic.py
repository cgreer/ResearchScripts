import bioLibCG
import updateDegPeaks as degPeaks
import sys
import blankIDs
from parRun import parRun
from checkExit import parClean


dRNAFN = sys.argv[1]
peakFN = sys.argv[2]
seqFN = sys.argv[3]

timer = bioLibCG.cgTimer()
timer.start()

#initialize oRNA database
print 'make blank IDs'
blankIDs.blankIDs(seqFN, dRNAFN)
print 'sequence'
degPeaks.updateSequence(dRNAFN, seqFN)
print 'tcc'
degPeaks.updateTcc(dRNAFN, peakFN)
print 'eLevel'
degPeaks.updateELevel(dRNAFN, '/home/chrisgre/smallLibs/siRNA/degradome/wigsk50')
print 'gSequence'
degPeaks.updateGSequence(dRNAFN)

print 'gScore'
parRun(50, 3, '/home/chrisgre/scripts/endoClip/updateDegPeaks.py', 'updateGScore', dRNAFN)
parClean(dRNAFN, 50)

