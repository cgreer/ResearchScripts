
from bioLibCG import *
import os
from optparse import OptionParser

def filterKnownMirs(gen = None, baseName = None):
	
	if gen is None:
		print 'Need assembly of genome!'
		return 1
	if baseName is None:
		print 'Need baseName of file family!'
		return 1
	
	#define databases to use (ensembl and mirBase for mouse, or for human)
	if gen == 'hg18':
		dbFileList = ['mirBaseHumanData/mirBaseData.dblColonDash', 'ensemblHumanData/ensemblData.dblColonDash']
		#mirBaseDataFileName = 'mirBaseHumanData/mirBaseData.dblColonDash'
		
	if gen == 'mm8':
		dbFileList = ['mirBaseMouseData/mirBaseData.dblColonDash', 'ensemblMouseData/ensemblData.dblColonDash']
		#mirBaseDataFileName = 'mirBaseMouseData/mirBaseData.dblColonDash'
	
	
	#This is the file where all the small RNA matches are...
	myDataFileName = '/u/home9/gxxiao/chrisgre/scripts/compareSmallRna/out/%s/%s.compared' % (baseName, baseName)
	#output file
	novelFileName = 'out/%s/%s.novel.mirs' % (baseName, baseName)

	chrs = []
	strands = ['1','-1']
	#define chr lists:
	if gen == 'hg18':
		for num in range(1,23):
			chrs.append(str(num))
		chrs.append('X')
		chrs.append('Y')

	elif gen == 'mm8':
		for num in range(1,20):
			chrs.append(str(num))
		chrs.append('X')
		chrs.append('Y')
	else:
		print 'Assembly not acceptable (currently only hg18 and mm8)'
		

	#create hitMap for mirBaseData
	print 'Creating hitMap for mirBase Data'

	baseHitMap = {} #format is hitMap['chr']['strand'] = {'Coordinate':'IDS'}

	for chr in chrs:
		baseHitMap[chr] = {}
		for strand in strands:
			baseHitMap[chr][strand] = {}
	
	for dbFile in dbFileList: #takes care of both ensembl and mirBase, I think overlapping dictionary entries are only counted once
			
		mirBaseDataFile = open(dbFile, 'r')

		for line in mirBaseDataFile:
			#strip data, generate list of coordinates for that mir, add coords to hitMap
			chromosome = line.strip().split(':')[0].split('chr')[1]
			strand = line.strip().split(':')[1]
			mirStart = line.strip().split(':')[2].split('-')[0]
			mirEnd = line.strip().split(':')[2].split('-')[1]
				
			baseID = 'NO_ID'
			#baseID = line.strip().split('\t')[8].split(' ')[1].split('"')[1]
			
			mirCoords = range(int(mirStart), int(mirEnd) + 1)
			
			for coord in mirCoords:
				if (strand == '+') and (chromosome in baseHitMap):
					baseHitMap[chromosome]['1'][coord] = baseID
				elif (strand == '-') and (chromosome in baseHitMap):
					baseHitMap[chromosome]['-1'][coord] = baseID
		mirBaseDataFile.close()


	#Compare my results with hitMap...
	print 'Comparing mirPipe results against mirBase hitmap'
	allIDs = []
	excludedIDs = []
	myDataFile = open(myDataFileName, 'r')

	for line in myDataFile:
		myMirID = '.'.join(line.strip().split(':')[1].split('.')[0:3])
		if myMirID not in allIDs: #then we'll check it
			allIDs.append(myMirID)
		else: #we've already checked this region
			continue
			
		myChromosome = line.strip().split(':')[1].split('.')[3].split('chr')[1]
		myMirStart = line.strip().split(':')[1].split('.')[2]
		myStrand = line.strip().split(':')[1].split('.')[4]
		
		#establish coordinates to check
		if myStrand == '1': #!!!CHECK FOR RIGHT NUMBER OF COORDS!!
			myMirCoords = range(int(myMirStart),int(myMirStart) + 22)
		else:
			myMirCoords = range(int(myMirStart) - 22, int(myMirStart))
		
		#print myMirCoords
		for myCoord in myMirCoords:
			if myCoord in baseHitMap[myChromosome][myStrand]:
				if myMirID not in excludedIDs:
					excludedIDs.append(myMirID)
					print 'took out %s' % myMirID
	myDataFile.close()

	#Now take out the IDS from all IDS and output to file/stats.
	totalIDs = len(allIDs)
	for ID in excludedIDs:
		allIDs.remove(ID)

	# Now that you have the IDs, go back and get the coordinates -> dblColonDash
	'''myDataFile = open(myDataFileName, 'r')

	for line in myDataFile:
		myMirID = '.'.join(line.strip().split(':')[1].split('.')[0:3])
		if myMirID not in allIDs: #then we'll check it
			allIDs.append(myMirID)
		else: #we've already checked this region
			continue
			
		myChromosome = line.strip().split(':')[1].split('.')[3].split('chr')[1]
		myFrameStart = line.strip().split(':')[1].split('.')[5]
		myStrand = line.strip().split(':')[1].split('.')[4]
		
		if myStrand = '1':
			start = int(myFrameStart)
			end = int(myFrameStart) + 110
		else:
			end = int(myFrameStart)
			start = end - 110
	'''
	#Make directory for files to go into
	filePath = 'out/%s' % baseName
	if baseName in os.listdir('out/'):
		#Delete all contents
		for file in os.listdir(filePath):
			os.remove('%s/%s' % (filePath,file) )
	else:
		os.mkdir(filePath)

	novelFile = open(novelFileName, 'w')

	allIDs.sort()
	for ID in allIDs:#!!! I want all the information here though, not just the id!
		novelFile.write('%s\n' % ID)

	novelFile.close()

	print '%s/%s IDs were NOVEL!' % (len(allIDs), totalIDs)

if __name__ == "__main__":
	import sys
	parser = OptionParser()
	parser.add_option("-b", "--baseName", dest="baseName", help="base name of file (LIThuman-s3k8b17) or other")
	parser.add_option("-g", "--genome", dest="genome", help="the genome assembly (e.g. hg18 for human)")

	(options, args) = parser.parse_args()
	
	filterKnownMirs(baseName = options.baseName, gen = options.genome)







		
