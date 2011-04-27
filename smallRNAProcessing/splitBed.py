#split bed file into corresponding chromosome files
import bioLibCG as cg

def splitBed(fN, organism):
	
	file = open(fN, 'r')
	file.readline() #skip first line
	
	#open chromosome files:
	chromFiles = {}
	if organism == 'human':
		for chrom in cg.humanChromosomes:
			chromFiles[chrom] = open(fN + '.' + chrom + '.wig', 'w')
			chromFiles[chrom].write('track type=bedGraph\n')
	elif organism == 'mouse':
		for chrom in cg.mouseChromosomes:
			chromFiles[chrom] = open(fN + '.' + chrom + '.wig', 'w')
			chromFiles[chrom].write('track type=bedGraph\n')
	else:
		print 'organism not supported'
		return 1
	
	for line in file:
		chr = line.split('\t')[0]
		if chr in chromFiles:
			chromFiles[chr].write(line)

if __name__ == "__main__":
	import sys
	
	splitBed(sys.argv[1], sys.argv[2])
		
