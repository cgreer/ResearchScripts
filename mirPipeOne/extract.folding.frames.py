#Given a tsv database with the conserved kmer among X number of species =>
#  output tsv database with sliding frames within a window sequence region.


#imports
import GenomeFetch
from optparse import OptionParser

#THESE SHOULD BE PASSED

def extractFrames(gnms, inF, outF, expL, winL, fraS):
	
	#defaults
	if not gnms: gnms = 'hg18;mm8;rn4;canFam2'
	if not inF:
		print 'Error: input file needed'
		return False
	if not outF:
		print 'Error: output file needed'
		return False
	if not expL: expL = 100
	if not winL: winL = 110
	if not fraS: fraS = 3
	
	genomes = gnms.split(';')
	inFile = inF
	outFile = outF
	expandLength = int(expL)
	windowLength = int(winL)
	frameStep = int(fraS)
	
	#initial variables
	kmersToExpand = [] # this will hold all the kmer sequences and their info
	kmerID = 0

	#extract the data from the tsv file

	lineFile = open(inFile, 'r')

	for line in lineFile:
		
		dataLine = line.strip().split('\t')
		
		#create kmer dict data entity for each genome
		kmerID = kmerID + 1
		for i,genome in enumerate(genomes):
			kmer = {}
			kmer['sequence'] = dataLine[0]
			kmer['chromosome'] = dataLine[1+i*3]
			kmer['strand'] = int(dataLine[2+i*3])
			kmer['start'] = int(dataLine[3+i*3])
			kmer['genome'] = genome
			kmer['ID'] = kmerID		

			#now add to expand list
			kmersToExpand.append(kmer)

	lineFile.close()


	################################################
	##create 100bp expanded sequence for each kmer##
	################################################
	#check the size
	if not kmersToExpand:
		print "No conserved Kmers in this region, quitting..."
		return 0
	else:
		print "....number of conserved kmers: %s" % int(len(kmersToExpand)/4)
	
	for kmer in kmersToExpand:
		if kmer['strand'] == 1:
			expandStart = kmer['start'] - expandLength #check for +1/-1 difference? !!!	
			expandEnd = kmer['start'] + len(kmer['sequence']) + expandLength -1  # why is the minus one here?
		else:
			expandStart = kmer['start'] - len(kmer['sequence']) + 1 - expandLength
			expandEnd = kmer['start'] + expandLength
		
		genomeFetch = GenomeFetch.GenomeFetch(kmer['genome'])
		kmer['window'] = genomeFetch.get_seq_from_to(kmer['chromosome'], expandStart, expandEnd, kmer['strand'])
	print '....window length: %s' % len(kmersToExpand[-1]['window'])

	#######################################
	##Create sliding frames for each kmer##
	#######################################
	
	for kmer in kmersToExpand:
		kmer['slidingFrames'] = []
		
		#set the initial frame end
		frameEnd = windowLength
		while frameEnd < len(kmer['window']):
			frameStart = frameEnd - windowLength
			kmer['slidingFrames'].append(kmer['window'][frameStart:frameEnd])

			frameEnd = frameEnd + frameStep

	print '....frame length: %s' % len(kmersToExpand[-1]['slidingFrames'][0])

	###############################
	##Make an RNAfold input file###
	###############################

	foldFile = open(outFile, 'w')

	for kmer in kmersToExpand:
		for i,frame in enumerate(kmer['slidingFrames']):
			#calculate framestart here !!! NOT ABOVE
			if kmer['strand'] == 1:
				frameStart = kmer['start'] - expandLength + i*frameStep
			else:
				frameStart = kmer['start'] + expandLength - i*frameStep
			
			#write out
			desc = '> %s.%s.%s.%s.%s.%s.%s.%s\n' % (kmer['ID'], i, kmer['genome'], kmer['chromosome'], kmer['strand'], kmer['sequence'], kmer['start'], frameStart)
			foldFile.write(desc)
			frameLine = '%s\n' % frame
			foldFile.write(frameLine)
	foldFile.write('@\n')

	foldFile.close()

if __name__ == "__main__":
	import sys
	parser = OptionParser()
	parser.add_option("-i", "--infile", dest="inFile", help="input file")
	parser.add_option("-o", "--outfile", dest="outFile", help="output file")
	parser.add_option("-g", "--genomes", dest="genomes", help="string of genome builds seperated by semicolon")
	parser.add_option("-w", "--wlength", dest="windowLength", help="Length of window to extract frames from")
	parser.add_option("-f", "--flength", dest="frameLength", help="Length of each frame")
	parser.add_option("-s", "--fstep", dest="frameStep", help="number of nt to skip between each frame")
	(options, args) = parser.parse_args()
	
	extractFrames(gnms = options.genomes, inF = options.inFile, outF = options.outFile, winL = options.windowLength, expL = options.frameLength, fraS = options.frameStep)
	
	
		


		
