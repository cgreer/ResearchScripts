import bioLibCG as cg
import cgConfig as c

def getTccFromSamLine(line):
	'''SAM has odd formatting at top'''
	try:
		lineSplit = line.strip().split('\t')
		chrom = lineSplit[2]
		strand = lineSplit[1]
		if strand == '16':
			strand = '-1'
		else:
			strand = '1'
		start = int(lineSplit[3])
		end = start + len(lineSplit[9])
		
		return cg.makeTcc(chrom, strand, start, end)
	except:
		print 'Warning: line failed parsing'
		print line.strip()
		return None
		
def getTccFromBowtieLine(line):
	chrom = line.strip().split('\t')[2]
	strand = line.strip().split('\t')[1]
	if strand == '+':
		strand = '1'
	else:
		strand = '-1'
	start = int(line.strip().split('\t')[3])
	end = start + len(line.strip().split('\t')[4])
	
	return cg.makeTcc(chrom, strand, start, end)

def getTccFromUCSCLine(line):
	'''format may change...'''
	chrom = line.strip().split('\t')[1]
	strand = line.strip().split('\t')[6]
	if strand == '+':
		strand = '1'
	else:
		strand = '-1'
	start = int(line.strip().split('\t')[2])
	end = int(line.strip().split('\t')[3])
	
	return cg.makeTcc(chrom, strand, start, end)
	
def getTccFromBedLine(line):
	chrom = line.strip().split('\t')[0]
	strand = line.strip().split('\t')[5]
	if strand == '+':
		strand = '1'
	else:
		strand = '-1'
	start = int(line.strip().split('\t')[1])
	end = int(line.strip().split('\t')[2])
	
	return cg.makeTcc(chrom, strand, start, end)

def returnParserFunction(option):
	'''interface'''
	if option == 'Bowtie':
		return getTccFromBowtieLine
	elif option == 'UCSC':
		return getTccFromUCSCLine
	elif option == 'Bed':
		return getTccFromBedLine
	elif option == 'Sam':
		return getTccFromSamLine

def writeWigFromHitDict(hitDict, assembly, name, directory = None):
	
	mConf = c.getConfig('Main.conf')
	if not directory: directory = mConf.conf['wigs']
	if not name: name = cg.getBaseFileName(name, naked = True)
	lDict = cg.returnChromLengthDict(assembly)
	
	cg.clearDirectory(directory, overwrite = False)
	#write results to wig file
	for chrom in hitDict:
		for strand in hitDict[chrom]:
			
			oF = open(directory + '/%s.%s.%s.wig' % (name, chrom, strand), 'w')
			oF.write('track type=bedGraph name=%s.%s.%s\n' % (name, chrom, strand))
			
			#print '  sorting'
			#print hitDict[chrom]
			chromEnd = lDict[chrom] #
			hitDict[chrom][strand][chromEnd] = 0
			keys = hitDict[chrom][strand].keys()
			keys.sort()
			
			#print '  writing blocks'
			prevVal = 0
			prevCoord = 0
			blockStart = 0
			blockEnd = 1
			for key in keys:
				val = hitDict[chrom][strand][key]
				
				if prevCoord == key - 1: 
					if val == prevVal:#should be combined
						blockEnd = key + 1
					else: #no zero block
						#write old block
						oF.write('%s\t%s\t%s\t%s\n' % (chrom, blockStart, blockEnd, prevVal)) #!make it a float value?
						#start new block
						blockStart = key
						blockEnd = key + 1
						
				else:
					#write old block
					oF.write('%s\t%s\t%s\t%s\n' % (chrom, blockStart, blockEnd, prevVal))
					#write zero block
					oF.write('%s\t%s\t%s\t%s\n' % (chrom, blockEnd, key, 0))
					#start new block
					blockStart = key
					blockEnd = key + 1
				
				prevVal = val
				prevCoord = key
			oF.close()

def makeWig(fN, assembly, format = None, name = None):
	
	'''format assumes bowtie
	suitible for medium mapped files.
	takes longer.'''
	#assume bowtie
	if not format: format = 'Bowtie'
	parserFunction = returnParserFunction(format)
	if not name: name = cg.getBaseFileName(fN, naked = True)
	lDict = cg.returnChromLengthDict(assembly)
	
	
	for chrom in lDict:
		if not chrom in cg.acceptableChroms: continue
		for strand in ['1', '-1']:
			f = open(fN, 'r')
			#create hitmap of chrom and strand
			print chrom, strand, 'hitmap'
			hitDict = {}
			for line in f:
				
				lChrom, lStrand, start, end = cg.tccSplit(parserFunction(line))
				lStrand = str(lStrand)
				start = int(start)
				end = int(end)
				if chrom == lChrom and strand == lStrand:
					for i in range(start, end + 1):
						try:
							hitDict[i] += 1
						except KeyError:
							hitDict[i] = 1
			
			#write results to wig file
			writeWigFromHitDict(hitDict, assembly)
			


def makeWigMem(fN, assembly, format = None, name = None, directory = None):
	'''format assumes bowtie
	suitible for small mapped files.'''
	
	if not name: name = cg.getBaseFileName(fN, naked = True)
	if not format: format = 'Bowtie'
	parserFunction = returnParserFunction(format)
	
	lDict = cg.returnChromLengthDict(assembly)
	f = open(fN, 'r')
	f.readline() #header...file might not have one but its one read...
	
	#create hitmap of chrom and strand
	hitDict = {} #format = chr: { strand : { coord : value 
	for line in f:
		try:
			lChrom, lStrand, start, end = cg.tccSplit(parserFunction(line))
		except AttributeError:
			continue
		lStrand = str(lStrand)
		start = int(start)
		end = int(end)
		if lChrom in cg.acceptableChroms:
			
			#wig for degradome
			if lStrand == '1':
				i = start + 20
			else:
				i = start
				
			try:
				hitDict[lChrom][lStrand][i] += 1
			except KeyError:
				if lChrom not in hitDict:
					hitDict[lChrom] = {}
				if lStrand not in hitDict[lChrom]:
					hitDict[lChrom][lStrand] = {}
				hitDict[lChrom][lStrand][i] = 1
			'''
			
			for i in range(start, end):
				try:
					hitDict[lChrom][lStrand][i] += 1
				except KeyError:
					if lChrom not in hitDict:
						hitDict[lChrom] = {}
					if lStrand not in hitDict[lChrom]:
						hitDict[lChrom][lStrand] = {}
					hitDict[lChrom][lStrand][i] = 1
			'''		
	f.close()
	
	#write results to wig file
	writeWigFromHitDict(hitDict, assembly, name, directory)

def mixWig(directory, assembly, name = None):
	'''Does it by chromosome --> faster, less memory'''
	
	if not name: name = 'Merge'
	#gather all chromosomes
	chromList = []
	for fN in cg.recurseDir(directory, end = '.wig'):
		chrom = cg.getBaseFileName(fN).strip().split('.')[-3]
		if chrom not in chromList:
			chromList.append(chrom)
	
	print chromList
	
	for chrom in chromList:
		
		print chrom
		#Gather all the values from all the files
		hitDict = {} # chrom : { strand : coord
		for fN in cg.recurseDir(directory, end = '.wig'):
			fChrom = cg.getBaseFileName(fN).strip().split('.')[-3]
			if fChrom != chrom: continue
			print  '  ', fN, fChrom
			f = open(fN, 'r')
			f.readline() #header
			strand = cg.getBaseFileName(fN).strip().split('.')[-2]
			for line in f:
				
				lChrom, start, end, val = (line.strip().split('\t'))
				start, end, val = int(start), int(end), int(val)
				if val < 1: continue
				#print start, end, val
				for i in range(start, end):
					try:
						hitDict[lChrom][strand][i] += val
					except (KeyError,TypeError):
						if not lChrom in hitDict:
							hitDict[lChrom] = {}
						if not strand in hitDict[lChrom]:
							hitDict[lChrom][strand] = {}
						hitDict[lChrom][strand][i] = val
		
		#write results to wig file
		writeWigFromHitDict(hitDict, assembly, name, directory)


if __name__ == "__main__":
	import sys
	
        cg.submitArgs(makeWigMem, sys.argv)
	#cg.submitArgs(mixWig, sys.argv)

				
				

	
