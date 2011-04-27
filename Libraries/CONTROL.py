import cgPlot
import cgSort
import compareData as compare
import stepVectorScan as sv
import cgPeaks

#cgPeaks test
'''
tcc = 'chr1:1:1102538:1103319'
stretch = cgPeaks.stretch(tcc, 'agoProfile.conf')
stretch.createContBlocks()

sortedKeys = stretch.profile.keys()
sortedKeys.sort()

for key in sortedKeys:
	print key, stretch.profile[key]

for block in stretch.blocks:
	print block
'''


'''
#histogram test
a = [5,5,5,5,5,5,4,4,4,3,3,2,1]
cgPlot.plotHistogram(a, name = 'testHistogram')
'''

#svCoord test

tcc = 'chr1:1:1103298:1103298'

eLevels = sv.svCoord([tcc], 'siPeaks.conf')
for coord in eLevels:
	print coord, eLevels[coord]

'''
tcc = 'chr1:1:1102539:1103319'
eLevels = sv.profileAroundPoint(tcc, 50, 'agoProfile.conf')
sorted = eLevels.keys()
sorted.sort()

for coord in sorted:
	print coord, eLevels[coord]
'''

#Plotting
'''
tccList = compare.tccFileToList('/home/chrisgre/dataSources/known/Human/tccs/mirBaseHuman.gff.tcc', 0)
for tcc in tccList:
	cgPlot.plotASProfile(tcc, 'H2.conf', directory = '/home/chrisgre/ASPlots', min = 5)
'''

#Sorting
'''
sl = [[7, 'micro', 3.0],[7, 'micro', 2.0], [3, 'ele', 7.0]]
pl = [2,0] #the end is most important
cl = [3,0]
ol = [False, False]
print cgSort.sortLines(sl, pl, cl, ol)
'''
