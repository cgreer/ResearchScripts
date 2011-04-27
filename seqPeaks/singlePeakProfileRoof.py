#given tcc, return best peak combo
import bioLibCG as cg
import cgConfig as c
import wigValue
import compareData as compare

#init
mConf = c.cgConfig('Main.conf')
conf = c.cgConfig()
pRange = 30

tccList = compare.tccFileToList('mouseKnownMirs.tcc', 0)
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

print 'profiling peaks'
stretchLengths = []
rNextValue = {}
for i in range(0,5):
	rNextValue[i] = []


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
			#print '  peak added', i
			tccPeaks.append(i)
	
	if not len(tccPeaks) > 0:
		continue
		
	#get highest peak of peaks, this will be the one that is profiled...
	hVal = 0
	hPeak = 0
	for i in tccPeaks:
		thisVal = peaks[chrom][strand][i]
		if thisVal > hVal:
			hVal = thisVal
			hPeak = i
	
	if not hVal > 30:
		print '  expression not high enough'
		continue
	
	#Now profile it +/- X bp left and right...
	cProfile = {}
	zeroPoint = hPeak
	zeroVal = peaks[chrom][strand][zeroPoint]
	errors = False
	for i in range(0,pRange):
		#get values, ratio for points --> update profile
		hCoord = zeroPoint + i
		try:
			hVal = wigValue.getWigValue(chrom, strand, hCoord)
		except:
			print 'Index Messed up', chrom, strand, hCoord
			errors = True
			break
		hRatio = float(hVal)/float(zeroVal)
		
		lCoord = zeroPoint - i
		try:
			lVal = wigValue.getWigValue(chrom, strand, lCoord)
		except:
			print 'Index Messed up', chrom, strand, lCoord
			errors = True
			break
		lRatio = float(lVal)/float(zeroVal)
		
		cProfile[i] = hRatio
		cProfile[-i] = lRatio
		
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
				if stretch > highest:
					highest = stretch
					endFinal = i - 1
					startFinal = startCurrent
					startCurrent = None
				else:
					startCurrent = None
			stretch = 0
	
	stretchLengths.append(highest)
	
	#get next values
	if highest > 14 and highest < 26:
		for i in range(0,5):
			if endFinal + i < 30:
				rNextValue[i].append(cProfile[endFinal + i])
			if startFinal - i > 0:
				rNextValue[i].append(cProfile[startFinal - i])
	

print timer.split()

outFile = open('mouseStretchLengths', 'w')
outFile.write('stretchLength\n')
for length in stretchLengths:
	outFile.write('%s\n' % length)
outFile.close()

for i in rNextValue:
	outFile = open('mouseStretchNexts.%s' % i, 'w')
	outFile.write('next\n')
	for val in rNextValue[i]:
		outFile.write('%s\n' % val)
	outFile.close()


