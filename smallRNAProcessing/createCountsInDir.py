import bioLibCG as cg
import cgRnaSeq

dir = '/home/chrisgre/smallLibs/zebrafish.embryo.GSE22068.20620952'


QNames = cg.recurseDir(dir, end = '.fastq') #fastqfile names
slNames = cg.recurseDir(dir, end = '.txt') #single sequence per line file names
fastNames = cg.recurseDir(dir, end = '.fa')
fastNames.extend(cg.recurseDir(dir, end = '.fna'))

for QFileName in QNames:
	print 'Creating counts for file', QFileName
	cgRnaSeq.createCountFileFastQ(QFileName)

for slName in slNames:
	print 'Creating counts for file', slName
	cgRnaSeq.createCountFileSL(slName)


for fastName in fastNames:
	print 'Creating counts for file', fastName
	cgRnaSeq.createCountFileFasta(fastName)
