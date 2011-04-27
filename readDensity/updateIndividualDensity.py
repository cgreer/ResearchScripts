####!!!I THINK THIS IS BUGGY, DOESN't make sense with chromosomes and strands!
import cgConfig
import bioLibCG as cg

mainConf = cgConfig.cgConfig('Main.conf')
conf = cgConfig.cgConfig()


wigFiles = cg.recurseDir(mainConf.conf['smallPath'], end = '.mapped.1.wig')
wigFiles.extend(cg.recurseDir(mainConf.conf['smallPath'], end = '.mapped.-1.wig'))
idDict = {}
for wigFileN in wigFiles:
	print wigFileN
	
	wigStrand = wigFileN.split('.')[?]

	for chrom in cg.mouseChromosomes:	
		print '  chrom', chrom	
		#init
		wigFile = open(wigFileN, 'r')
		mirFile = open(conf.conf['results'], 'r')
	
		#get rid of header
		wigFile.readline()
	
		print '  populating hitmap'
		#populate hitmap
		wigMap = {}
		for line in wigFile:
			wigChrom = line.strip().split('\t')[0]
			value = int(line.strip().split('\t')[3].split('.')[0])
			if value > 0 and wigChrom == chrom:
			
				start = int(line.strip().split('\t')[1])
				end = int(line.strip().split('\t')[2])
				for i in range(start, end):
					wigMap[i] = value
		wigFile.close()
	
		print '  calculating hits for mature seqs'
		#calculate total hits per mature
		inList = []
		for line in mirFile:
			mTcc = line.strip().split('\t')[1]
			mirID = line.strip().split('\t')[0].split('.')[0]
			if mirID in inList: continue
		
			chrom = mTcc.split(':')[0]

			if wigChrom != chrom: continue

			totalHits = 0
			for i in range(int(mTcc.split(':')[2]), int(mTcc.split(':')[3]) + 1): 
				if i in wigMap:
					totalHits += wigMap[i]
	
			#update
			if mirID in idDict:
				if wigFileN in idDict[mirID]:
					idDict[mirID][wigFileN] += totalHits
				else:
					idDict[mirID][wigFileN] = totalHits
			else:#initialize the dictionary, initialize value
				idDict[mirID] = {}
				idDict[mirID][wigFileN] = totalHits
	
		mirFile.close()

print 'Writing New File'
#write new results file
outFile = open('individual.densities.data', 'w')
for mirID in idDict:
	outFile.write( mirID + '\n')
	for lib in idDict[mirID]:
		outFile.write('\t%s\t%s\n' % (lib, idDict[mirID][lib]))

outFile.close()


		
			
				
		
	
	
				
