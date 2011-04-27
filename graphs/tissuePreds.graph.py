#how much density for tissue?
import cgConfig
import bioLibCG as cg

mainConf = cgConfig.cgConfig('Main.conf')
conf = cgConfig.cgConfig()

predFile = open(conf.conf['resultsSorted'], 'r') #must be sorted!!!

mirIDs = [x.strip().split('\t')[0].split('.')[0] for x in predFile]

#make lib:tissue dict
metaDict = cg.getMetaFileDict(mainConf.conf['metaFileName'])
#lib is getBaseFileName(key, naked = True): tisue is metaDict[key][3]  


#make lib:hits dict
densityFile = open('/home/chrisgre/scripts/readDensity/individual.densities.data', 'r')
tissueHits = {}
cID = 'NONE'
for line in densityFile:
	if line.startswith('\t'): #lib: hit
		l = line.strip().split('\t')[0]
		hits = int(line.strip().split('\t')[1])
		tissueHits[cID][l] = hits 
	else:
		cID = line.strip()
		tissueHits[cID] = {}

tissueHist = {}#tissue: hits
for mirID in mirIDs:
	if mirID in tissueHits:
		for smallLib in tissueHits[mirID]:
			smallName = cg.getBaseFileName(smallLib, naked = True)
			if smallName in metaDict:
				if metaDict[smallName][1] == 'mouse':
					if len(metaDict[smallName]) > 3:
						t = metaDict[smallName][3]
						if t in tissueHist:
							tissueHist[t] += tissueHits[mirID][smallLib]
						else:
							tissueHist[t] = tissueHits[mirID][smallLib]

print tissueHist
