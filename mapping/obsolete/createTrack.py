import HTSeq
import bioLibCG as cg
import cgConfig as c
import subprocess, time, os
import clusterCheck
import cgSort

def updateWigLength(fN, assembly):
	
	chromLengths = cg.returnChromLengthDict(assembly)
	f = open(fN, 'r')
	header = f.readline() #header
	
	lineDict = {} # chr : []
	for line in f:
		chrom = line.split('\t')[0]
		if chrom not in lineDict:
			lineDict[chrom] = []
			lineDict[chrom].append(line)
		else:
			lineDict[chrom].append(line)
	f.close()
	
	for chrom in lineDict:
		print 'extending', chrom
		
		chromLength = chromLengths[chrom]
		lastChromLine = lineDict[chrom][-1]
		lastValue = int(lastChromLine.split('\t')[2])
		print lastValue
		if lastValue < chromLength:
			print '  ', lastValue, chromLength
			lineDict[chrom].append('%s\t%s\t%s\t0.000000\n' % (chrom, lastValue, chromLength))
			print '  updated'
		
	f = open(fN, 'w')
	f.write(header)
	for chrom in lineDict:
		f.writelines(lineDict[chrom])
	f.close()

def createTrack(fName, organism):

	if organism == 'human':
		chroms = cg.humanChromosomes
		assembly = 'hg19'
	elif organism == 'mouse':
		chroms = cg.mouseChromosomes
		assembly = 'mm9'
	elif organism == 'zebrafish':
		chroms = cg.zebrafishChromosomes
		assembly = 'danRer6'

	alignment_file = HTSeq.BowtieReader(fName)
	cvg = HTSeq.GenomicArray(chroms, stranded=True, typecode='i')
	for alngt in alignment_file:
		if alngt.aligned:
			cvg.add_value( 1, alngt.iv ) #iv is the genomic interval..

	bedNamePos = fName + '.1.' + 'wig'
	bedNameNeg = fName + '.-1.' + 'wig'

	cvg.write_bedgraph_file(bedNamePos, "+" )
	cvg.write_bedgraph_file(bedNameNeg, "-" )
	
	#Now extend it and sort it.
	updateWigLength(bedNamePos, assembly)
	updateWigLength(bedNameNeg, assembly)
	
	#Now Sort it.
	cgSort.wigSort(bedNamePos)
	cgSort.wigSort(bedNameNeg)

def createMTrack(dirName):
	'''merge all mapped tracks in directory and create a single wig file'''
	
	fileList = cg.recurseDir(dirName, end = '.out')
	
	chroms = cg.humanChromosomes
	print 'Making Bed File vectors'
	cvg = HTSeq.GenomicArray(chroms, stranded=True, typecode='i')
	for fName in fileList:
		print fName
		alignment_file = HTSeq.BowtieReader(fName)
		for alngt in alignment_file:
			if alngt.aligned:
				try:
					cvg.add_value( 1, alngt.iv ) #iv is the genomic interval..
				except KeyError:
					pass

	bedNamePos = dirName + '/Merge.' + 'hg19' + '.1.wig'
	bedNameNeg = dirName + '/Merge.' + 'hg19' + '.-1.wig'
	
	print 'Writing Bed File'
	cvg.write_bedgraph_file(bedNamePos, "+" )
	cvg.write_bedgraph_file(bedNameNeg, "-" )

	#Now extend it
	updateWigLength(bedNamePos, 'hg19')
	updateWigLength(bedNameNeg, 'hg19')
	
	#Now Sort it.
	cgSort.wigSort(bedNamePos)
	cgSort.wigSort(bedNameNeg)

def createMultiTrack(dirName, organism):
	'''merge all mapped tracks in directory and create a single wig file'''
	mainConf = c.cgConfig('Main.conf')
	metaFileName = mainConf.conf['metaFileName']
	
	fileList = []
	for file in cg.recurseDir(dirName, end = '.mapped'):
						
		#check if mouse or human SHOULD PUT INTO A STD FUNCTION FOR META FILE
		#check if mouse or human
		baseFName = cg.getBaseFileName(file, naked= True)
		
		metaDict = cg.getMetaFileDict(metaFileName)
		
		org = 'None'
		if baseFName in metaDict:
			if metaDict[baseFName][1] == 'NONE':
				print '  NO ORG KNOWN FOR', file
				continue
			elif not metaDict[baseFName][1] == organism:
				print '  NOT ORGANISM RUNNING', file
				continue
			else:
				org = metaDict[baseFName][1]
				print '  USING ORG', org, file
			
		#check if there is an organism, must check due to files not in metaFile
		if org == 'None':
			print '  NO org (not in meta file)', file
			continue
		
		#only make wig file for organism asked for
		if not org == organism:
			continue
		
		#if it is right organism and has mapped file then add
		fileList.append(file)
	
	
	#make merged wig
	if organism == 'human':
		chroms = cg.humanChromosomes
		assembly = 'hg19'
	elif organism == 'mouse':
		chroms = cg.mouseChromosomes
		assembly = 'mm9'
	elif organism == 'zebrafish':
		chroms = cg.zebrafishChromosomes
		assembly = 'danRer6'
	
	print 'Making Bed File vectors'
	cvg = HTSeq.GenomicArray(chroms, stranded=True, typecode='i')
	for fName in fileList:
		alignment_file = HTSeq.BowtieReader(fName)
		for alngt in alignment_file:
			if alngt.aligned:
				cvg.add_value( 1, alngt.iv ) #iv is the genomic interval..

	bedNamePos = dirName + '/Merge.' + organism + '.1.wig'
	bedNameNeg = dirName + '/Merge.' + organism + '.-1.wig'
	
	print 'Writing Bed File'
	cvg.write_bedgraph_file(bedNamePos, "+" )
	cvg.write_bedgraph_file(bedNameNeg, "-" )

	#Now extend it
	updateWigLength(bedNamePos, assembly)
	updateWigLength(bedNameNeg, assembly)
	
	#Now Sort it.
	cgSort.wigSort(bedNamePos)
	cgSort.wigSort(bedNameNeg)

def createMultiTrackDir(dirName, organism):
	'''THIS DIFFERS FROM ABOVE BECAUSE IT DOESN't REQUIRE META INFO
	IT JUST MAKES A MERGED WIG FOR EVERYTHING IN THE DIRECTORY'''
	mainConf = c.cgConfig('Main.conf')
	
	fileList = []
	for file in cg.recurseDir(dirName, end = '.mapped'):
		fileList.append(file)
				
	#make merged wig
	if organism == 'human':
		chroms = cg.humanChromosomes
		assembly = 'hg19'
	elif organism == 'mouse':
		chroms = cg.mouseChromosomes
		assembly = 'mm9'
	elif organism == 'zebrafish':
		chroms = cg.zebrafishChromosomes
		assembly = 'danRer6'
	
	print 'Making Bed File vectors'
	cvg = HTSeq.GenomicArray(chroms, stranded=True, typecode='i')
	for fName in fileList:
		alignment_file = HTSeq.BowtieReader(fName)
		for alngt in alignment_file:
			if alngt.aligned:
				cvg.add_value( 1, alngt.iv ) #iv is the genomic interval..

	bedNamePos = dirName + '/Merge.' + organism + '.1.wig'
	bedNameNeg = dirName + '/Merge.' + organism + '.-1.wig'
	
	print 'Writing Bed File'
	cvg.write_bedgraph_file(bedNamePos, "+" )
	cvg.write_bedgraph_file(bedNameNeg, "-" )

	#Now extend it
	updateWigLength(bedNamePos, assembly)
	updateWigLength(bedNameNeg, assembly)
	
	#Now Sort it.
	cgSort.wigSort(bedNamePos)
	cgSort.wigSort(bedNameNeg)


def createTrackInDir(dirName):
	'''Every Q function has a corresponding shell script
	Make wig file for all mapped files, for all organisms'''
	
	wrapperShell = '/home/chrisgre/scripts/mapping/createTrack.sh'
	
	mainConf = c.cgConfig('Main.conf')
	metaFileName = mainConf.conf['metaFileName']

	for file in cg.recurseDir(dirName, end = '.mapped'):
						
		#check if mouse or human
		baseFName = cg.getBaseFileName(file)
		baseFName = baseFName.split('.')[0]
		
		metaDict = cg.getMetaFileDict(metaFileName)
		
		org = 'None'
		if baseFName in metaDict:
			if metaDict[baseFName][1] == 'NONE':
				print '  NO ORG KNOWN FOR', file
				continue
			else:
				org = metaDict[baseFName][1]
				print '  USING ORG', org, file
				
		#check if there is an organism, must check due to files not in metaFile
		if org == 'None':
			print '  NO org (not in meta file)', file
			continue
			
		while True:
			#submit job if there are less than ten
			if clusterCheck.queryNumJobsQ('chrisgre') < 1000:
				#subprocess.Popen(['qsub', '-q', 'xiao', '-l', 'mem=4G', '-V', '-o', mainConf.conf['outLog'], '-e', mainConf.conf['errorLog'], wrapperShell, file, org ])
				subprocess.Popen(['qsub', '-V', '-o', mainConf.conf['outLog'], '-e', mainConf.conf['errorLog'], wrapperShell, file, org ])
				#time.sleep(.5) #give it time to update qstat
				break
			else:#wait 10 secs...
				time.sleep(20)

if __name__ == "__main__":
	import sys
	
	cg.submitArgs(createMTrack, sys.argv)
	#createTrack(sys.argv[1], sys.argv[2])

