import bioLibCG
import siRnaPredict as si
import sys
import updateContigs
import updateDuplicatesMultiTcc as udmt
import cgAlignment
import updateMismatchAndMiddle as mm
import filtering

#this is the directory where the alignment file has been placed
runName = sys.argv[1]

'''Pipeline assumes that full alignment of small to target seqs has been done, but not truncated'''
'''Also assumes that the results files have been made...'''
oRNADir = runName + '/oRNA'
aDir = runName + '/aDir'


timer = bioLibCG.cgTimer()
timer.start()

print 'filtering targets'
filtering.filterTargets(oRNADir, aDir, True, 1, 1, .50)
print timer.split()

