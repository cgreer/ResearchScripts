#given tcc, return best peak combo
import bioLibCG as cg
import cgConfig as c
import wigValue
import compareData as compare

#init
mConf = c.cgConfig('Main.conf')
conf = c.cgConfig()

tccList = ['chr3:-1:96042576:96042685']
#tccList = compare.tccFileToList('mouseKnownMirs.tcc', 0)
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
for tcc in tccList:
	print tcc
	tccPeaks = []
	chrom = cg.ss(tcc, ':')[0]
	strand = cg.ss(tcc, ':')[1]
	start = int(cg.ss(tcc, ':')[2])
	end = int(cg.ss(tcc, ':')[3])
	
	#get all peaks
	for i in range(start, end + 1):
		if i in peaks[chrom][strand]:
			print '  peak added', i
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
				
				#calc distance
				dist = y - x
				
				#calc ratio
				ratio = float(xval)/float(yval)
				if ratio < float(1):
					ratio = float(1)/ratio
				
				#height of biggest peakFiles in pair
				if xval > yval:
					max = xval
				else:
					max = yval
				
				#get hairpin ratio (point in between two peaks)
				#might want to only do this if distance is within zone...
				hPoint = x + dist/2
				hVal = wigValue.getWigValue(chrom, strand, hPoint)
				if hVal == 0:
					hRatio = max
				else:
					hRatio = float(max)/(hVal)
				
				peakCombos.append([tcc,x,y,dist,ratio,max,hRatio])
				#print '  ', peakCombos[-1]
	
	#find best combo...
	topCombo = None
	for i, combo in enumerate(peakCombos):
		dist = int(combo[3])
		hRatio = int(combo[6])
		if (22 < dist < 68) and (hRatio > 2): #distance (22,57) and hRatio(5) are 2 sds!
			if not topCombo:
				topCombo = combo
			else: #take one with highest hRatio
				if combo[6] > topCombo[6]:
					print 'switched', topCombo, combo
					topCombo = combo
	
	if topCombo:
		bestCombos.append(topCombo)
		print bestCombos[-1]

print timer.split()

'''
#output results to a file for R
outFile = open('mousePeaks.R.data', 'w')
outFile.write('tcc\tpeakOne\tpeakTwo\tdistance\tratio\tmax\thRatio\n')
for peak in bestCombos:
	outFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n' % (peak[0], peak[1],peak[2],peak[3],peak[4],peak[5],peak[6]))
outFile.close()
'''
