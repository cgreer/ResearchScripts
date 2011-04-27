import bioLibCG as cg

fileName = '/home/chrisgre/smallLibs/blelloch.mouse.ES.pmid.18923076/GSM314553-10948.counts'

file = open(fileName, 'r')
outFile = open(fileName + '.fasta', 'w')

for line in file:
	seq = line.strip().split('\t')[0]
	if '.' in seq: #clean check
		continue
		
	count = int(line.strip().split('\t')[1])
	i = 0
	if i < count:
		outFile.write('> \n')
		outFile.write(seq + '\n')
		i = i + 1

