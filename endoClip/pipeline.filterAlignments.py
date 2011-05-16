import bioLibCG
import sys
from parRun import parRun
from checkExit import parClean
import subprocess

#fN and Init
aFN = sys.argv[1] 
aFilteredFN = sys.argv[2]

timer = bioLibCG.cgTimer()
timer.start()

#filter targets
print 'filtering targets'
parRun(30, 3, '/home/chrisgre/scripts/endoClip/filteringFlat.py', 'filterTargetsInPlace', aFN, 'True', '1', '1', '.50')
parClean(aFN, 30)
print timer.split()

#make db smaller
print 'truncating db'
subprocess.Popen(['/home/chrisgre/scripts/endoClip/truncate.filterAlignments.sh', aFN, aFilteredFN]).wait()
print timer.split()




