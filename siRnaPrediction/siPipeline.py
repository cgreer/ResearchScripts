import siRnaPredict as si
import cleanup
import subprocess

peakFName =     'small.peak.data'
smallUpdateFN = 'small.degradome'
cleanFN =       'small.degradome.small'
exonUpdateFN =  'small.degradome.small.clean'
mismatchFN =    'small.degradome.small.clean.degOverlap'
clean2FN =      'small.degradome.small.clean.degOverlap.mismatch'
sortFN =        'small.degradome.small.clean.degOverlap.mismatch.clean'

smallDirectory = '/home/chrisgre/smallLibs/siRNA/small/sortedChroms'

print 'Getting highest level, offset, middle level'
si.predictSiRNA(peakFName, 'siDegradome.conf')
print 'Getting small expression level'
si.updateSmallExpression(smallUpdateFN, 'siPeaks.conf')
print 'cleaning'
cleanup.cleanup(cleanFN)
print 'Getting transcript level'
si.exonOverlap(exonUpdateFN)
print 'Getting mismatch levels'
si.markMismatches(mismatchFN, smallDirectory)
print 'Cleaning up again, and sorting'
cleanup.cleanup2(clean2FN)
subprocess.Popen(['./sortIt.sh', sortFN], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
