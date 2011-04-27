

def updateLength(inFile, outFile):

	f = open(inFile, 'r')
	fOut = open(outFile, 'w')
	for line in f:
		ls = line.strip().split('\t')
		length = int(ls[2])
		newLength = length + 1
		
		fOut.write('%s\t%s\t%s\n' % (ls[0], ls[1], newLength))

	f.close()
	fOut.close()

if __name__ == '__main__':
	import sys
	updateLength(sys.argv[1], sys.argv[2])
