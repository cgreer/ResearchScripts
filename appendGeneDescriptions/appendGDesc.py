import getHairpins
import cgGenes
import cgConfig as c
import bioLibCG as cg

mConf = c.getConfig('Main.conf')
geneSetFolder = mConf.conf['geneSetsHuman']

fN = '/home/chrisgre/projects/NoncodingHuman/results/NChuman-s3k8b17.results.sorted.introns.sorted'
cHairs = getHairpins.getHairpins(fN)

ensGenes = cgGenes.createGeneSetFromFile(geneSetFolder + '/ensemblAllTranscripts.tsv')

cDesc = {} #CID:gDesc
for CID in cHairs:
	tcc = cHairs[CID]
	
	cDesc[CID] = "NONE"
	
	overlappingGenes = ensGenes.geneOverlaps([tcc])
	if len(overlappingGenes) > 0:
		print overlappingGenes[0].type
		cDesc[CID] = overlappingGenes[0].type

f = open(fN, 'r')
newLines = []
for line in f:
	CID = line.strip().split('\t')[7]
	newLines.append(cg.appendToLine(line, cDesc[CID], 16))
f.close()

f = open(fN + '.FINAL', 'w')
f.writelines(newLines)
f.close()
		
		
	
	
		
	
