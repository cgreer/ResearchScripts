#Convert Sequence/Count file to psuedo fastq so it can be parsed by my fastq script.
#Some of the read data is in the format of a sequence and the count (# times that sequence was read)
##I need to convert it to fastq so I can parse it and have the smallRNA comparison work well.

baseName = 'undiff_mapped.txt'
libFileName = '/u/home9/gxxiao/apps/projects/small.rna.libs/marra.human.ES/%s' % baseName
libFile = open(libFileName, 'r')

outFileName = '%s.fastq' % baseName
outFile = open(outFileName, 'w')

for line in libFile:
	
	#get rid of header
	if line.startswith('#') or line.startswith('ID') or line.startswith('Se'):
		print 'skipped [%s]' % line
		continue
	else:
		sequence = line.strip().split('\t')[1]
		outFile.write('@ %s\n%s\n' % (baseName, sequence))

