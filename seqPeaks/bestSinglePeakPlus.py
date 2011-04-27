#given tcc, return best peak combo
import bioLibCG as cg
import cgConfig as c
import wigValue
import compareData as compare
import getHairpins
import math
import stepVectorScan
import cgPeaks

def findPeaks(pType, cName = None):
	
	#init
	mConf = c.cgConfig('Main.conf')
	conf = c.getConfig(cName)

	if pType == 'E':
		predName = conf.conf['resultsExonsSorted']
	else:
		predName = conf.conf['resultsIntronsSorted']
	
	print predName
	#make CID:hairpin:peak dictionary
	cHairs = getHairpins.getHairpins(predName)
	peakDict = {}
	for CID in cHairs:
		peakDict[CID] = [cHairs[CID],'None']
		

	timer = cg.cgTimer()
	timer.start()

	#put peaks in memory
	print 'Creating peak data'
	peaks = {} # chr:peak:value
	for CID in cHairs:
		chrom, strand, start, end = cg.tccSplit(cHairs[CID])
		tcc = cHairs[CID]
		
		#init dictionary
		if chrom not in peaks:
			peaks[chrom] = {}
		
		if strand not in peaks[chrom]:
			peaks[chrom][strand] = {}
		
		#create peaks for tcc and add to peak dictionary
		stretch = cgPeaks.stretch(tcc, cName)
		stretch.createPeaks()
		for peakCoord in stretch.peaks:
			peaks[chrom][strand][peakCoord] = 0
	print timer.split()

	print 'finding best combos'
	bestCombos = []
	aPass = 0
	bPass = 0
	cPass = 0
	numT = 0
	for CID in peakDict:
		cgFlag = False
		if CID == '538':cgFlag = True
		
		tcc = peakDict[CID][0]
		#print tcc
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
				
								
				#scan a 30 bp range around this point and find the best roof...
				pRange = 30
				rTcc = cg.makeTcc(chrom, strand, x, x + 1)
				
				#quickly get max value...kinda a long way to do it but whatever
				cProfile = stepVectorScan.profileAroundPoint(rTcc, 1, cName, ratio = False)
				xval = cProfile[0]
				max = xval
				highestValueCoord = x
				
				#now make profile for roof...
				cProfile = stepVectorScan.profileAroundPoint(rTcc, pRange, cName, ratio = True)
				
				
				
				#now get highest stretch length and the rNext coord.
				minVal = .80
				highest = 0
				stretch = 0
				startCurrent = None
				startFinal = None
				endFinal = None
				for i in range(1 - pRange, pRange):
					if cProfile[i] > minVal:
						stretch += 1
						if startCurrent == None:
							startCurrent = i
					else:
						if stretch > 0:
							if stretch > highest: #stretch ended and was higher than previous
								highest = stretch
								endFinal = i - 1
								startFinal = startCurrent
								startCurrent = None
							else:
								startCurrent = None
						stretch = 0
				
				#get +/- 4 value...
				val = [1.0, 1.0]
				if (startFinal) and (endFinal):
					low = startFinal - 4
					high = endFinal + 4
					if low > (1 - pRange):
						if high < pRange:
							val[0] = float(cProfile[startFinal - 4])
							val[1] = float(cProfile[endFinal + 4])
				
				#fill in other details...
				y = 'S'
				dist = 'S'
				ratio = 'S'
				
				peakCombos.append([tcc,x,y,dist,ratio,max,highest,val])
				#print '  ', peakCombos[-1]
		
		#find best combo...
		topCombo = None
		for combo in peakCombos:
			roofLength = combo[6]
			dropValue = combo[7][0]
			if combo[7][1] > dropValue:
				dropValue = combo[7][1]
			
			#print roofLength, dropValue
			if 14 < roofLength < 26:
				if 0.0 < dropValue < 0.2:
					#pick one with rooflength nearest 20:
					if topCombo:
						if (math.fabs(22 - roofLength)) < (math.fabs(22 - topCombo[6])):
							topCombo = combo
					else:
						topCombo = combo
		
		if topCombo:
			peakDict[CID][1] = topCombo
			bestCombos.append(topCombo)
			print bestCombos[-1]
		else:
			#print 'None'
			pass

	print timer.split()


	#now update predFile (SLOT 13)
	predFile = open(predName, 'r')
	newLines = []
	for line in predFile:
		CID = cg.ss(line)[7]
		if peakDict[CID][1] == 'None':
			peakInfo = 'None'
		else:
			peakInfo = '%s:%s:%s:%s:%s:%s' % (str(peakDict[CID][1][1])[-3:], 'S', str(peakDict[CID][1][4]).split('.')[0], peakDict[CID][1][5],peakDict[CID][1][6], peakDict[CID][1][7])
		newLines.append(cg.appendToLine(line, peakInfo, 13))
	predFile.close()

	predFile = open(predName, 'w')
	predFile.writelines(newLines)
	predFile.close()

if __name__ == "__main__":
	import sys

	if len(sys.argv) > 2:
		findPeaks(sys.argv[1], sys.argv[2])
	else:
		findPeaks(sys.argv[1])	
