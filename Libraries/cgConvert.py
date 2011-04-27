import bioLibCG as cg

#order lists that tell how the data is formatted
#chrom, strand, start, end, strandType, sep

bedOrder = [0,5,1,2,1,'\t' , 'chr']
UCSCOrder = [1,6,2,3,1,'\t', 'chr']
gffOrder = [0,6,3,4,1,'\t', '']
tccOrder = [0,1,2,3,0,':', 'chr']

def returnOrderList(format):
	if format == 'Bed':
		return bedOrder
	elif format == 'UCSC':
		return UCSCOrder
	elif format == 'Gff':
		return gffOrder
	elif format == 'Tcc':
		return tccOrder
	else:
		print 'Format given does not exist'
		return 1
		

def convertFile(fN, inFormat, outFormat, oFN = None):
	
	#get input order and extract info
	iO = returnOrderList(inFormat)
	oO = returnOrderList(outFormat)
        
	
	f = open(fN, 'r')
	newLines = []
	for line in f:
		ls = line.strip().split(iO[5])
		chrom, strand, start, end = ls[iO[0]], ls[iO[1]], ls[iO[2]], ls[iO[3]]
	        
                #switch to appropriate chromosome type if needed
                if len(chrom) == 1:
                        chrom = oO[6] + chrom

		#switch strand if need be
		if oO[4] == 0:
			if strand == '1' or strand == '+':
				strand = '1'
			else:
				strand = '-1'
		else:
			if strand == '1' or strand == '+':
				strand = '+'
			else:
				strand = '-'
		
		#construct new Line
		newLine = '\n'
		newLine = cg.appendToLine(newLine, chrom, oO[0], sep = oO[5])
		newLine = cg.appendToLine(newLine, strand, oO[1], sep = oO[5])
		newLine = cg.appendToLine(newLine, start, oO[2], sep = oO[5])
		newLine = cg.appendToLine(newLine, end, oO[3], sep = oO[5])
		
		newLines.append(newLine)
	f.close()
	
	#output file
	f = open(fN + '.' + outFormat, 'w')
	f.writelines(newLines)
	f.close()

if __name__ == "__main__":
	import sys
	cg.submitArgs(convertFile, sys.argv)
