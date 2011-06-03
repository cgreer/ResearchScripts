import bioLibCG as cg
from bioLibCG import ss
import cgConfig as c
import cgIndex
#init



'''Given Tcc List
Feach TCC, goto correct line in index
scan values till end of tcc
goto next tcc'''


def scanVectorsHist(tccList, cName):
	'''Given tcc list --> scan wig files and get histogram values
	can be modified to do single/total values...
	THIS USES INDEXES!!! = BAD...'''
	
	conf = c.getConfig(cName)
	org = conf.conf['organism']
	mConf = c.getConfig('Main.conf')
	wigDir = mConf.conf['wig%s' % org]

	
	timer = cg.cgTimer()
	timer.start()
	histDict = {} # tcc: [list values]
	for tcc in tccList:
		theSplit = ss(tcc, ':')
		chrom, strand, tccStart, tccEnd = theSplit[0], theSplit[1],int(theSplit[2]),int(theSplit[3])
		
		#goto correct fild, correct line in index
		
		fN = wigDir + '/Merge.%s.%s.wig.%s.wig' %  (org.lower(),strand,chrom)
		fNindex = wigDir + '/Merge.%s.%s.wig.%s.wig.index' % (org.lower(),strand,chrom)
		
		#print timer.split()
		#get line in index file
		iFile = open(fNindex, 'r')
		startByte = 'None'
		for line in iFile:
			beg = int(cg.ss(line)[1])
			end = int(cg.ss(line)[2])
			
			if beg <= tccStart < end:
				startByte = int(cg.ss(line)[0]) 
				#print 'INDEX', line.strip()
				break
		iFile.close()
		
		#print timer.split()
		#grab value
		f = open(fN, 'r')
		f.seek(startByte, 0)
		
		stop = False
		for line in f:
			#print 'Line:', line.strip()
			lBeg = int(cg.ss(line)[1])
			lEnd = int(cg.ss(line)[2])
			lValue = int(cg.ss(line)[3].split('.')[0])
			
			if tccStart > lBeg:
				lBeg = tccStart
			if tccEnd < lEnd:
				lEnd = tccEnd
				stop = True
			#print timer.split()

			for i in range(lBeg, lEnd):
				try:
					histDict[tcc].append(lValue)
				except KeyError: #just for zero...so you don't have to if every time...
					histDict[tcc] = [lValue]
			if stop: break
	
		f.close()
		#print timer.split()
	return histDict

def scanVectorsSingleCoord(tccList, cName):
	'''Given tcc list --> scan wig files and coord:value...
	'''
	
	conf = c.getConfig(cName)
	org = conf.conf['organism']
	mConf = c.getConfig('Main.conf')
	wigDir = mConf.conf['wig%s' % org]

	timer = cg.cgTimer()
	timer.start()
	coordDict = {} # tcc: [list values]
	for tcc in tccList:
		theSplit = ss(tcc, ':')
		chrom, strand, tccStart, tccEnd = theSplit[0], theSplit[1],int(theSplit[2]),int(theSplit[3])
		
		#goto correct fild, correct line in index
		
		fN = wigDir + '/Merge.%s.%s.wig.%s.wig' %  (org.lower(),strand,chrom)
		fNindex = wigDir + '/Merge.%s.%s.wig.%s.wig.index' % (org.lower(),strand,chrom)
		
		#print timer.split()
		#get line in index file
		iFile = open(fNindex, 'r')
		startByte = 'None'
		for line in iFile:
			beg = int(cg.ss(line)[1])
			end = int(cg.ss(line)[2])
			
			if beg <= tccStart < end:
				startByte = int(cg.ss(line)[0]) 
				#print 'INDEX', line.strip()
				break
		iFile.close()
		
		#print timer.split()
		#grab value
		f = open(fN, 'r')
		f.seek(startByte, 0)
		
		stop = False
		for line in f:
			#print 'Line:', line.strip()
			lBeg = int(cg.ss(line)[1])
			lEnd = int(cg.ss(line)[2])
			lValue = int(cg.ss(line)[3].split('.')[0])
			
			if tccStart > lBeg:
				lBeg = tccStart
			if tccEnd < lEnd:
				lEnd = tccEnd
				stop = True
			#print timer.split()

			for i in range(lBeg, lEnd):
				coordDict[i] = lValue
				
			if stop: break
	
		f.close()
	return coordDict

def profileAroundPoint(zeroPoint, span, cName, ratio = False, ratioCoord = None):
	'''span is +/- that number...  if you put in 30, you'll get 
	-29 to 29 around 0
	zeroPoint must be in tcc format with point at start
	ratio will return all points around zero point as a ratio of zero point
	so zeroPoint = 1.0 and the rest will be fractions of that...'''
	chrom, strand, zPoint, end = cg.tccSplit(zeroPoint)
	rStart = zPoint - span
	rEnd = zPoint + span
	
	rTcc = cg.makeTcc(chrom, strand, rStart, rEnd)
	scanDict = scanVectorsSingleCoord([rTcc], cName)
	
	#reorient so that zero is at zero 
	#find position of zero
	returnDict = {}
	if not ratio:
		for i in range(1-span, span):
			returnDict[i] = scanDict[zPoint + i]
	else:
		if ratioCoord:
			zeroVal = scanDict[ratioCoord]
		else:
			zeroVal = scanDict[zPoint]
		if zeroVal == 0: zeroVal = 1 
		for i in range(1-span, span):
			r = float(scanDict[zPoint + i])/float(zeroVal)
			returnDict[i] = r
	
	return returnDict

def scanVectorsOrganism(tccList, config = None):
	'''Given tcc list --> scan Organism wig files and coord:value...
	'''
	
	config = c.getConfig(config)
	
	coordDict = {} # tcc: [list values]
	for tcc in tccList:
		chrom, strand, tccStart, tccEnd = cg.tccSplit(tcc)
		
		#print 'Checking Tcc'	
		org = config.conf['organism']
		mConf = c.getConfig('Main.conf')
		wigDir = mConf.conf['wig%s' % org]
		fN = wigDir + '/Merge.%s.%s.wig.%s.wig' %  (org.lower(),strand,chrom)	
		#print 'Checking Index'
		#goto correct line in index
		fIndex = cgIndex.lineIndex(fN, header = True) #!!!there actually is a header...have to deal with this...
		fIndex.passCheckFunction(cgIndex.wigCheckFunction)
		fIndex.binarySearch(tcc) #places file pointer at beginning of tcc as beginning
		
		stop = False
		for line in fIndex.file:
			
			#print 'Line:', line.strip()
			lBeg = int(cg.ss(line)[1])
			lEnd = int(cg.ss(line)[2])
			lValue = int(cg.ss(line)[3].split('.')[0])
			
			if tccStart > lBeg:
				lBeg = tccStart
			if tccEnd < lEnd:
				lEnd = tccEnd
				stop = True
			#print timer.split()

			for i in range(lBeg, lEnd):
				coordDict[i] = lValue
				
			if stop: break
	return coordDict

def scanVectorsFile(fN, tccList):
	'''Given tcc list --> scan wig files and return coord:value...
	'''	
	timer = cg.cgTimer()
	timer.start()
	coordDict = {} # tcc: [list values]
	for tcc in tccList:
		chrom, strand, tccStart, tccEnd = cg.tccSplit(tcc)
		
		#goto correct line in index
		fIndex = cgIndex.lineIndex(fN, header = True) #!!!there actually is a header...have to deal with this...
		fIndex.passCheckFunction(cgIndex.wigCheckFunction)
		fIndex.binarySearch(tcc) #places file pointer at beginning of tcc as beginning
				
		stop = False
		for line in fIndex.file:
			#print 'Line:', line.strip()
			lBeg = int(cg.ss(line)[1])
			lEnd = int(cg.ss(line)[2])
			lValue = int(cg.ss(line)[3].split('.')[0])
			
			if tccStart > lBeg:
				lBeg = tccStart
			if tccEnd < lEnd:
				lEnd = tccEnd
				stop = True
			#print timer.split()

			for i in range(lBeg, lEnd):
				coordDict[i] = lValue
				
			if stop: break
	
		#fIndex.close()
	return coordDict

def svCoord(tccList, config = None):
	'''Given tcc list --> scan Organism wig files and coord:value...
	'''
	
	#init
	config = c.getConfig(config)
	org = config.conf['organism']
	wigDir = config.conf['wigSetDir']
	wigSetName = config.conf['wigSetName']
	splitIntoChroms = config.conf['wigChromSplit']
	if splitIntoChroms == 'True':
		splitIntoChroms = True
	else:
		splitIntoChroms = False

	coordDict = {} # tcc: [list values]
	for tcc in tccList:
		chrom, strand, tccStart, tccEnd = cg.tccSplit(tcc)
		
		if splitIntoChroms:
			fN = wigDir + '/%s.%s.%s.wig' %  (wigSetName, chrom, strand)
		else:
			fN = wigDir + '/Merge.%s.%s.wig' % (org.lower(), strand)
		
		fIndex = cgIndex.lineIndex(fN, header = True)
		fIndex.passCheckFunction(cgIndex.wigCheckFunction)
		fIndex.binarySearch(tcc) #places file pointer at beginning of tcc as beginning
		
		stop = False
		for line in fIndex.file:
			
			#print 'Line:', line.strip()
			lBeg = int(cg.ss(line)[1]) + 1
                        #print 'lBeg', lBeg
			lEnd = int(cg.ss(line)[2])
                        #print 'lEnd', lEnd
                        #print '--'
			lValue = int(cg.ss(line)[3].split('.')[0])
			
			if tccStart > lBeg:
				lBeg = tccStart
			if tccEnd < lEnd:
				lEnd = tccEnd
				stop = True
			#print timer.split()

			for i in range(lBeg, lEnd + 1):
				coordDict[i] = lValue
				
			if stop: break
		fIndex.close() #close the file and the index after use...

	return coordDict
        
def profileAroundPoint(zeroPoint, span, cName, ratio = False, ratioCoord = None):
	'''span is +/- that number...  if you put in 30, you'll get 
	-29 to 29 around 0
	zeroPoint must be in tcc format with point at start
	ratio will return all points around zero point as a ratio of zero point
	so zeroPoint = 1.0 and the rest will be fractions of that...'''
	chrom, strand, zPoint, end = cg.tccSplit(zeroPoint)
	rStart = zPoint - span
	rEnd = zPoint + span
	
	rTcc = cg.makeTcc(chrom, strand, rStart, rEnd)
	scanDict = svCoord([rTcc], cName)
	
	#reorient so that zero is at zero 
	#find position of zero
	returnDict = {}
	if not ratio:
		for i in range(1-span, span):
			returnDict[i] = scanDict[zPoint + i]
	else:
		if ratioCoord:
			zeroVal = scanDict[ratioCoord]
		else:
			zeroVal = scanDict[zPoint]
		if zeroVal == 0: zeroVal = 1 
		for i in range(1-span, span):
			r = float(scanDict[zPoint + i])/float(zeroVal)
			returnDict[i] = r
	
	return returnDict
