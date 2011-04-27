#Turn fastq into fasta

fN = ''

def qToFasta(fN):
	
	#go to first line with @, 
	f = open(fN, 'r')
	while not f.readline() == '@':
		pass
	
	outF = open(fN + '.fasta', 'w')
	i = 0
	for line in f:
		if i%4 == 0:
			outF.write('>\n%s\n' % line.strip())

		i = i + 1

if __name__ == "__main__":
	import sys
	
	qToFasta(sys.argv[1])
