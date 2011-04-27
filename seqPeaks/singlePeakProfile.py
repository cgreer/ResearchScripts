#given tcc, return best peak combo
import bioLibCG as cg
import cgConfig as c
import wigValue
import compareData as compare

#init
mConf = c.cgConfig('Main.conf')
conf = c.cgConfig()
pRange = 100



tccList = ['chr3:-1:96042576:96042685', 'chr3:-1:96042576:96042685']
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

#build profile data structure
profile = {}
for i in range(0,pRange):
	profile[i] = []
	profile[-i] = []

print 'profiling peaks'
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
	zeroPoint = hPeak
	zeroVal = peaks[chrom][strand][zeroPoint]
	errors = False
	for i in range(1,pRange):
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
		
		profile[i].append(hRatio)
		profile[-i].append(lRatio)
	
	if not errors:
		profile[0].append(1.0) #get max number of items to tally...


print timer.split()

list = profile.keys()
list.sort()
'''
for i in list:
	print i, profile[i]
'''

#output results to a file for R
outFile = open('peakProfile.data', 'w')
newList = [str(x) for x in list]
outFile.write('\t'.join(newList) + '\n')
newLines = []
print len(profile[0])
for j in range(0, len(profile[0])):
	#print 'lineNum', j
	newLine = []
	for i in range(1 - pRange,pRange):
		newLine.append(str(profile[i][j]))
	newLine.append('\n')
	newLines.append('\t'.join(newLine))
	#print newLines[-1]
	
outFile.writelines(newLines)
outFile.close()

