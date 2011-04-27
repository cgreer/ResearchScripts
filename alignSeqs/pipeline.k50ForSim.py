import bioLibCG
import siRnaPredict as si
import sys
import updateContigs
import updateDuplicatesMultiTcc as udmt
import cgAlignment
import updateMismatchAndMiddle as mm

#this is the directory where the alignment file has been placed
runName = sys.argv[1]

'''Pipeline assumes that full alignment of small to target seqs has been done, but not truncated'''
'''Also assumes that the results files have been made...'''
oRNADir = runName + '/oRNATest'
aDir = runName + '/aDir'

alignTrunFN = runName + '/all.aligned.0.5.filtered'
seqFN = '/home/chrisgre/scripts/alignSeqs/oRNA.sequences'
tFN = '/home/chrisgre/scripts/alignSeqs/dRNA.results.updated'

timer = bioLibCG.cgTimer()
timer.start()

#initialize oRNA database
print 'initializing oRNA'
si.updateID(oRNADir, seqFN)
print timer.split()

#initializing alignments
print 'Initializing alignments'
print '...loading alignments'
cgAlignment.loadAlignments2(aDir, alignTrunFN)
print timer.split()
print '...appending T Info'
cgAlignment.appendTInfo(aDir, tFN)
print '...appending Tran Info'
cgAlignment.appendTranInfo(aDir, tFN)
print timer.split()


print '...updating paired interactions: centered mismatches and center expression'
mm.markCenterExpression(aDir, 'siDegradome.conf')
mm.markMismatchedPairs(aDir)
print timer.split()

print 'Linking oRNA to Targets'
#update the initial targets for small rnas based off of only local alignment
si.updateTargetIDs(oRNADir, aDir)
print timer.split()

print 'DONE with', runName
