import cgConfig as c
import bioLibCG as cg
import cgIndex

def scanSequence(seqList, dirName):
	'''Given list of sequences --> get all reads that have sequence
	'''
	
	fileNames = cg.recurseDir(dirName, end = '.sequence')
	if len(fileNames) > 1:
		print fileNames
		print 'there is more than one sequence file in this directory'
		return 1
	else:
		fN = fileNames[0]
	
	#for seq in seqList:
	seq = seqList
	fIndex = cgIndex.lineIndex(fN, header = False)
	fIndex.passCheckFunction(cgIndex.mapSequenceCheckFunction)
	fIndex.binarySearch(seq) #places file pointer at beginning of sequence line
	
	#extend and report
	fIndex.extendUp(seq)
	finalReads = []
	for line in fIndex.file:
		if fIndex.checkFunction(seq, line) == 0:
			finalReads.append(line.strip())
		else:
			return finalReads


def scanCoord(tcc, dirName):
	
	fileNames = cg.recurseDir(dirName, end = '.starts')
	
	
	#get name of file for index
	chrom, strand, start, end = cg.tccSplit(tcc)
	nameCheck = '%s.%s' % (chrom, strand)
	fN = 'None'
	for fileName in fileNames:
		if nameCheck in fileName: fN = fileName
	if fN == 'None': 
		print 'No Index file for', nameCheck
		return 0
        
	
	fIndex = cgIndex.lineIndex(fN, header = False)
	fIndex.passCheckFunction(cgIndex.mapStartCheckFunction)
	fIndex.binarySearch(tcc, skipEnd = True) #places file pointer at beginning of sequence line
        
        #Check if you need to move down one line
        checkLine = fIndex.getLineFromByte(fIndex.currentByte)
        fIndex.passCheckFunction(cgIndex.mapStartRangeCheckFunction) #Note i'm passing now, but it is also used in extending
        if fIndex.checkFunction(tcc, checkLine) != 0:
                fIndex.file.readline()
        fIndex.currentByte = fIndex.file.tell()

        #Now extend up until in range, down until in range --> return reads.
        fIndex.extendUp(tcc)
        
	finalReads = []
	for line in fIndex.file:
		if fIndex.checkFunction(tcc, line) == 0:
			finalReads.append(line.strip())
                else:
			return finalReads
        


if __name__ == "__main__":
	import sys
        for read in cg.submitArgs(scanCoord, sys.argv):
                print read
	
