#split preliminary predictions into exons and introns.
import bioLibCG as cg
from bioLibCG import ss
import cgConfig as c
import getHairpins
import compareData as compare

def splitExonsIntrons(cName = None):
	mConf = c.cgConfig('Main.conf')
	conf = c.getConfig(cName)
	
	#init
	organism = conf.conf['organism']
	minOverlap = 50
	cHairs = getHairpins.getHairpins() #CID: HAIRPIN
	exonList = compare.tccFileToList('%sExons.tcc' % organism, 0)
	hairpins = []
	for CID in cHairs:
		hairpins.append(cHairs[CID])
	
	print 'checking overlaps'
	#check which hairpins overlap exons and by how much
	exonOverlapped = compare.compareTwoTcc(hairpins, exonList, 1, amount = True)
	print '  ', len(exonOverlapped)
	
	print 'removing partial introns'
	#remove the ones that didn't overlap more than X:
	remList = []
	for tcc, oAmount in exonOverlapped:
		if oAmount < minOverlap:
			remList.append([tcc, oAmount])
	
	for item in remList:
		exonOverlapped.remove(item)
	print '  ', len(exonOverlapped), 'out of', len(cHairs.keys())
		
	#get CIDs of exons
	exonCIDs = []
	for tcc, oAmount in exonOverlapped:
		for CID in cHairs:
			if cHairs[CID] == tcc:
				exonCIDs.append(str(CID))
	
	
	#Open sorted predictions and write lines with CIDs to respective files
	predFile = open(conf.conf['resultsSorted'], 'r')
	exonFile = open(conf.conf['resultsSorted'] + '.exons', 'w')
	intronFile = open(conf.conf['resultsSorted'] + '.introns', 'w')
	for line in predFile:
		if line.split('\t')[7] in exonCIDs:
			exonFile.write(line)
		else:
			intronFile.write(line)
	predFile.close()
	exonFile.close()
	intronFile.close()

if __name__ == "__main__":
	import sys
	if len(sys.argv) > 1:
		splitExonsIntrons(sys.argv[1])
	else:
		splitExonsIntrons()
