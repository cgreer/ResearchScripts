#given tcc, return best peak combo
import bioLibCG as cg
import cgConfig as c
import wigValue
import compareData as compare
import getHairpins
import math

#init
mConf = c.cgConfig('Main.conf')
conf = c.cgConfig()
predName = conf.conf['resultsExonsSorted']

#make CID:hairpin:peak dictionary
cHairs = getHairpins.getHairpins(predName)
peakDict = {}
for CID in cHairs:
	peakDict[CID] = [cHairs[CID],'None']
	

timer = cg.cgTimer()
timer.start()

#put peaks in memory
print 'loading peak data'
peakFilesNames = cg.recurseDir(mConf.conf['wigMouse'], end = '.peaks')
peaks = {} # chr:peak:value
for pN in peakFilesNames:
	chrom = pN.strip().split('.')[4]
	strand = pN.strip().split('.')[2]
	
	#init dictionary
	if chrom not in peaks:
		peaks[chrom] = {}
	
	if strand not in peaks[chrom]:
		peaks[chrom][strand] = {}
	
	#get peaks and values and put in dictionary
	pFile = open(pN, 'r')
	for line in pFile:
		peaks[chrom][strand][int(line.strip().split('\t')[0])] = int(line.strip().split('\t')[1].split('.')[0])
print timer.split()

print 'finding best combos'
bestCombos = []
aPass = 0
bPass = 0
cPass = 0
numT = 0
for CID in peakDict:
	
	tcc = peakDict[CID][0]
	print tcc
	tccPeaks = []
	chrom = cg.ss(tcc, ':')[0]
	strand = cg.ss(tcc, ':')[1]
	start = int(cg.ss(tcc, ':')[2])
	end = int(cg.ss(tcc, ':')[3])
	
	#get all peaks
	for i in range(start, end + 1):
		if i in peaks[chrom][strand]:
			#print '  peak added', i
			tccPeaks.append(i)
	
	#Calculate parameters...
	pairStrings = [] #used to check if pair already added
	peakCombos = []
	for x in tccPeaks:
			
			#peak values
			xval = peaks[chrom][strand][x]
			max = xval
			highestValueCoord = x
			
			#get 6 points ratios
			samplePoints = []
			for diff in [7,15,20]:
				for mult in [1, -1]:
					newDiff = diff*mult
					val = wigValue.getWigValue(chrom, strand, highestValueCoord + newDiff)
					r = float(val)/float(max)
					samplePoints.append(r)
			
			y = 0
			dist = 0
			ratio = 0
			
			
			peakCombos.append([tcc,x,y,dist,ratio,max,samplePoints])
			#print '  ', peakCombos[-1]
	
	#find best combo...
	topCombo = None
	for combo in peakCombos:
		numT +=1
		sPoints = combo[6]
		max = combo[5]
		aFlag = False
		bFlag = False
		cFlag = False
		if sPoints[0] > .97 or sPoints[1] > .97: #flat points
			print 'passed A'
			aFlag = True
			aPass += 1
		if sPoints[2] < .1 or sPoints[3] < .1: #drop points
			print 'passed B'
			bFlag = True
			bPass += 1
		if sPoints[4] < .01 or sPoints[5] < .01: #zero points
			print 'passed C'
			cFlag = True
			cPass += 1
		
		#did it pass?
		if aFlag and bFlag and cFlag:
			if not topCombo:
				topCombo = combo
			else: #take one with highest max peak.
				if max > topCombo[5]:
					topCombo = combo
	
	if topCombo:
		peakDict[CID][1] = topCombo
		bestCombos.append(topCombo)
		print bestCombos[-1]
	else:
		print 'None'

print timer.split()
print 'Total', numT
print 'A', aPass
print 'B', bPass
print 'C', cPass


#output results to a file for R
outFile = open('mousePeaksResults.R.data', 'w')
outFile.write('tcc\tpeakOne\tpeakTwo\tdistance\tratio\tmax\thRatio\n')
for peak in bestCombos:
	outFile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (peak[0], peak[1],peak[2],peak[3],peak[4],peak[5]))
outFile.close()


#now update predFile (SLOT 13)
predFile = open(predName, 'r')
newLines = []
for line in predFile:
	CID = cg.ss(line)[7]
	if peakDict[CID][1] == 'None':
		peakInfo = 'None'
	else:
		peakInfo = '%s:%s:%s:%s' % (str(peakDict[CID][1][1])[-3:], str(peakDict[CID][1][2])[-3:], str(peakDict[CID][1][4]).split('.')[0], peakDict[CID][1][5])
	newLines.append(cg.appendToLine(line, peakInfo, 13))
predFile.close()

predFile = open(predName, 'w')
predFile.writelines(newLines)
predFile.close()
	
