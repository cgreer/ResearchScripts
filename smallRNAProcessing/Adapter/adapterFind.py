from optparse import OptionParser

def adapterFind(fName, oName, org = None):
	
	#init:
	if org is None:
		print 'setting org to human'
		org = 'human'
		
	#mir sequences: let-7, mir21
	if org == 'human':
		mirSeqs = ['TGAGGTAGTAGGTTGTATAGTT', 'TAGCTTATCAGACTGATGTTGA', 'CTTTTTGCGGTCTGGGCTTGC', 'CAAAGTGCTTACAGTGCAGGTAG']
	elif org == 'mouse':
		mirSeqs = ['TGAGGTAGTAGGTTGTATAGTT', 'TAGCTTATCAGACTGATGTTGA']
	
	for mirSeq in mirSeqs:
		reads = []
		append = reads.append
		qFile = open(fName, 'r')
		
		#Get lines with mir seq
		for line in qFile:
			if mirSeq in line:
				append(line.strip())
		qFile.close()
		
		#Count lines with each adapter sequence
		seqDict = {}
		for read in reads:
			if read.replace(mirSeq, ' ') in seqDict:
				seqDict[read.replace(mirSeq, ' ')] = seqDict[read.replace(mirSeq, ' ')] + 1
			else:
				seqDict[read.replace(mirSeq, ' ')] = 1
				
		#Output Distribution to outFile or std out:
		sortDict = {}
		if seqDict:
			for key in seqDict:
				sortDict[seqDict[key]] = key
		
		sortList = []
		for key in sortDict: #Now the counts...
			sortList.append(key)
		sortList.sort()
		sortList.reverse() #by highest
		
		#populate final list
		if len(sortList) > 10:
			r = 9
		else:
			r = len(sortList)
		
		finalList = []
		i = 0
		while i < r:
			finalList.append('%s\t%s' % (sortDict[sortList[i]], sortList[i] ))
			i = i + 1
		

		if oName is None:
			print '\nTop Adapters for', mirSeq
			for line in finalList:
				print line
		else:
			oFile = open(oName, 'r')
			oFile.write('\nTop Adapters for %s\n' % mirSeq)
			for line in finalList:
				oFile.write('%s\n' % line)
			oFile.close()

if __name__ == "__main__":
	import sys
	parser = OptionParser()
	parser.add_option("-i", dest="fName", help="name of file")
	parser.add_option("-o", dest="oName", help="name of output file")
	parser.add_option("-g", dest="org", help="organism (human or mouse)")

	(options, args) = parser.parse_args()
	
	adapterFind(fName = options.fName, oName = options.oName, org = options.org)


