

#imports
import GenomeFetch
import re
from optparse import OptionParser
from doubleloop import *

def UtoT(string):
	letters = list(string)
	for i,letter in enumerate(letters):
		if letter == 'U':
			letters[i] = 'T'
	return ''.join(letters)

def filterFrames(inF = None, mirLength = None, genomes = None, mirBasePairs = None, outF = None):
	
	#defaults
	if not inF:
		print 'Error: input file needed for filtering frames'
		return False
	if not outF:
		print 'Error: output file needed for filtering frames'
		return False

	if not mirLength: mirLength = 22
	if not genomes: genomes = 'hg18;mm8;rn4;canFam2'
	if not mirBasePairs: mirBasePairs = 16
	
	#cast correctly
	gnms = genomes.split(';')
	firstGenome = gnms[0]
	mirLength = int(mirLength)
	mirBasePairs = int(mirBasePairs)
	
	#open files
	print '....opening file'
	lineFile = open(inF, 'r')
		
	
	#######################
	#Make 22mers (or xmers)
	#######################
	print '....extracting potential mir sequences'	
	
	pidList = [] #this is used to do Genome inclusive checks.
	gMirs = {}
	for genome in gnms:
		gMirs[genome] = {} #will contain ids linked to lists with mirs in them.
	
	while True:
		
		frame = {} #make new frame
		#grab lines
		lineOne = lineFile.readline()
		if lineOne == '': break
		lineTwo = lineFile.readline()
		if lineTwo == '': break
		lineThree = lineFile.readline()
		if lineThree == '': break
		
		#parse three lines
		desc = lineOne.strip().split(' ')[1].split('.')
		frame['kmerID'] = desc[0]
		frame['frameID'] = desc[1]
		frame['genome'] = desc[2]
		frame['chromosome'] = desc[3]
		frame['strand'] = desc[4]
		frame['seed'] = desc[5].upper()
		frame['kmerStart'] = desc[6]
		frame['frameStart'] = desc[7]
		frame['sequence'] = UtoT(lineTwo.strip()) #second line
		frame['energy'] = float(lineThree[-8:-2]) #this might not be the best way to do it... !!!
		frame['structure'] = lineThree.split(' ')[0]
		
		
		#Filter single frame
		passEnergy = False
		passLoop = False
		MLCheck = True
		if frame['energy'] < -25.00:
			passEnergy = True
		
		loopMatch = re.search('[(][.]{2,100}[)]', frame['structure']) 	
		if loopMatch is not None: #means there was a match
			passLoop = True
		
		#double or more loop check
		try:
			if multiLoopCheck(frame['seed'], frame['sequence'], frame['structure']) == False:
				MLCheck = False
		except:
			print 'multiloop messed up'
			MLCheck = False
		
		if passEnergy and passLoop and MLCheck:
			mirBegin = 0
			mirEnd = mirLength
			while mirEnd <= len(frame['sequence']):
				mir = {}
				#sequence
				mir['sequence']  = frame['sequence'][mirBegin:mirEnd]
				mir['structure'] = frame['structure'][mirBegin:mirEnd]
				mirBegin = mirBegin + 1
				mirEnd = mirEnd + 1
				
				#Check if seed region is at 5' end
				seedLength = len(frame['seed'])
				if mir['sequence'][0:seedLength] == frame['seed'].upper():
					pass
				else: 
					continue
				
				#Check if contains at least 16 bp (or however many)
				if (mir['structure'].count(')') >= mirBasePairs) or (mir['structure'].count('(') >= mirBasePairs):
					pass
				else:
					continue
				
				#Check if has no loop inside of it ( '(' followed by ')' )
				firstLeft = mir['structure'].find('(')
				if firstLeft >= 0:
					if mir['structure'].find(')',firstLeft) >= 0:
						continue
	
				#If it passed all the tests, add it to mirs.
				mirID = frame['kmerID'] + '.' + frame['frameID']
				gen = frame['genome']
				mir['kmerID'] = frame['kmerID']
				mir['genome'] = frame['genome']
				mir['frameID'] = frame['frameID']
				mir['chromosome'] = frame['chromosome']
				mir['strand'] = frame['strand']
				mir['seed'] = frame['seed']
				mir['kmerStart'] = frame['kmerStart']
				mir['frameStart'] = frame['frameStart']
				mir['frameEnergy'] = frame['energy']
				#STILL MISSING WHERE THE FRAME STARTED CAN PROBABLY CALCULATE BASED OFF OF FRAME ID AND STEP... !!!

				#ADD
				gMirs[gen][mirID] = mir
				# add this prospective id to list
				if mirID not in pidList:
					pidList.append(mirID)
	lineFile.close()
	
	#Do cross genome check
	print '......amount of mirs before cross genome check: %s' % len(pidList)
	#Now collect all mirs whose tests passed for all 4 genomes
	
	finalMirs = {}
	for pid in pidList:
		flag = True
		
		#check if the mir with the id passed in all four genomes.
		for gen in gMirs.keys():
			if not gMirs[gen].has_key(pid):
				flag = False
		#if it passed, add it to final list.
		if flag:
			#finalMirs[pid] = []
			finalMirs[pid] = {}
			for gen in gMirs.keys():
				#finalMirs[pid].append(gMirs[gen][pid])
				finalMirs[pid][gen] = gMirs[gen][pid]

	print '...... # of prospective pre-mirs: %s' % len(finalMirs)
	
	outFile = open(outF, 'w')
	
	for id in finalMirs.keys():
		
		myStrand = finalMirs[id][firstGenome]['strand']
		myChromosome = finalMirs[id][firstGenome]['chromosome']
		myEnergy = finalMirs[id][firstGenome]['frameEnergy']
		myMirSeq = finalMirs[id][firstGenome]['sequence']
		
		#make standardized hairpin coords (triple colon and start always < end)
		if myStrand == '1':
			myFrameStart = finalMirs[id][firstGenome]['frameStart']
			myFrameEnd = int(myFrameStart) + 109 ###!!! a little lazy - technically this could change...but we haven't used anything besides 110 nt...
		else:
			myFrameEnd = finalMirs[id][firstGenome]['frameStart']
			myFrameStart = int(myFrameEnd) - 110 + 1 ###!!! check for +/- 1 error
		
		frameCoords = '%s:%s:%s:%s' % (myChromosome, myStrand, myFrameStart, myFrameEnd)
		
		#make standardized mature microRNA coords
		if myStrand == '1':
			myMirStart = finalMirs[id][firstGenome]['kmerStart']
			myMirEnd = int(myMirStart) + int(mirLength) -1 ###!!! a little lazy - technically this could change...but we haven't used anything besides 110 nt...
		else:
			myMirEnd = finalMirs[id][firstGenome]['kmerStart']
			myMirStart = int(myMirEnd) - int(mirLength) + 1 ###!!! check for +/- 1 error
		
		mirCoords = '%s:%s:%s:%s' % (myChromosome, myStrand, myMirStart, myMirEnd)
		
		desc = '%s.%s\t%s\t%s\t%s\t%s\n' % (id, myFrameStart, mirCoords, frameCoords, myMirSeq, myEnergy)
		outFile.write(desc)
	
	outFile.close()
	
	
if __name__ == "__main__":
	import sys
	parser = OptionParser()
	parser.add_option("-i", "--infile", dest="inFile", help="input file")
	parser.add_option("-o", "--outfile", dest="outFile", help="output file")
	parser.add_option("-g", "--genomes", dest="genomes", help="string of genome builds seperated by semicolon")
	parser.add_option("-m", "--mlength", dest="mirLength", help="Length of potential mirs")
	parser.add_option("-b", "--mirbp", dest="mirBasePairs", help="minimum # of bps in mirna (inversely, how little bulges)")
	(options, args) = parser.parse_args()
	
	filterFrames(inF = options.inFile, mirLength = options.mirLength, genomes = options.genomes, mirBasePairs = options.mirBasePairs, outF = options.outFile)
	
		


		
