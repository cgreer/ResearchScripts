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
	
	#get all combos
	pairStrings = [] #used to check if pair already added
	peakCombos = []
	for x in tccPeaks:
		for y in tccPeaks:
			if x != y: #don't make pair with self
				#check to see if we didn't already make this pair
				s1 = '%s:%s' % (x,y)
				s2 = '%s:%s' % (y,x)
				if (s1 in pairStrings) or (s2 in pairStrings):
					continue
				else:
					pairStrings.extend([s1,s2])
				
				#make sure x is lower of two...
				if x > y:
					x,y = y,x
				
				#peak values
				xval = peaks[chrom][strand][x]
				yval = peaks[chrom][strand][y]
				
				#highestPeak
				
				#calc distance
				dist = y - x
				
				#calc ratio
				ratio = float(xval)/float(yval)
				if ratio < float(1):
					ratio = float(1)/ratio
				
				#height of biggest peakFiles in pair
				if xval > yval:
					max = xval
					highestValueCoord = x
				else:
					max = yval
					highestValueCoord = y
				
				#get 7 points
				samplePoints = []
				samplePoints.append(wigValue.getWigValue(chrom, strand, highestValueCoord + 7))
				samplePoints.append(wigValue.getWigValue(chrom, strand, highestValueCoord - 7))
				samplePoints.append(wigValue.getWigValue(chrom, strand, highestValueCoord + 15))
				samplePoints.append(wigValue.getWigValue(chrom, strand, highestValueCoord - 15))
				samplePoints.append(wigValue.getWigValue(chrom, strand, highestValueCoord + 20))
				samplePoints.append(wigValue.getWigValue(chrom, strand, highestValueCoord - 20))
				'''
				#get hairpin ratio (point in between two peaks)
				#might want to only do this if distance is within zone...
				hPoint = x + dist/2
				hVal = wigValue.getWigValue(chrom, strand, hPoint)
				if hVal == 0:
					hRatio = max
				else:
					hRatio = float(max)/(hVal)
				'''
				
				peakCombos.append([tcc,x,y,dist,ratio,max,samplePoints])
				#print '  ', peakCombos[-1]
	
	#find best combo...
	topCombo = None
	for combo in peakCombos:
		dist = int(combo[3])
		sPoints = combo[6]
		if sPoints[0] > .97 and sPoints[1] > .97: #flat points
			print 'passed A'
			if sPoints[2] < .5 and sPoints[3] < .5: #drop points
				print 'passed B'
				if sPoints[4] < .2 and sPoints[5] < .2: #zero points
					print 'passed C'
					if (22 < dist < 57): #distance (22,57) and hRatio(5) are 2 sds!
						print 'passed distance'
						if not topCombo:
							topCombo = combo
						else: #take one with best peak distance seperation.
							if (math.fabs(37 - dist)) < (math.fabsabs(37 - topCombo[3])):
								print 'switched', topCombo, combo
								topCombo = combo
	
	if topCombo:
		peakDict[CID][1] = topCombo
		bestCombos.append(topCombo)
		print bestCombos[-1]
	else:
		print 'None'

print timer.split()


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
	
