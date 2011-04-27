import cgGenes
import compareData as compare
import cgConfig as c

cName = 'mm9.conf'
mConf = c.getConfig('Main.conf')
conf = c.getConfig(cName)
organism = conf.conf['organism']
geneSetFolder = mConf.conf['geneSets%s' % organism]
genes = cgGenes.createGeneSetFromFile(geneSetFolder + '/allTransciptsType.tsv')
peakTccs = compare.tccFileToList('peakData.500.mm9', 0)


tOverlaps = genes.transcriptOverlaps(peakTccs)
typeDict = {}
for transcript in tOverlaps:
	if transcript.type not in typeDict:
		typeDict[transcript.type] = 1
	else:
		typeDict[transcript.type] += 1

#count the amounts of each type for each transcript
amount = {}
for gene in genes.genes:
	for t in gene.transcripts:
		if t.type in amount:
			amount[t.type] += 1
		else:
			amount[t.type] = 1

print 'Total Peaks:', len(peakTccs)
print 'Total Peaks with Annotations', len(tOverlaps)
print 'Annotation Counts:'
for type in typeDict:
	percentage = float(typeDict[type])/float(amount[type])
	print ' ', type, typeDict[type], '\t', amount[type], str(percentage)
	


