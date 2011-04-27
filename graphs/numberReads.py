#generate number of reads per prediction.

import bioLibCG as cg
import cgConfig

conf = cgConfig.cgConfig()
resultsFile = open(conf.conf['resultsSorted'], 'r')
resultsFile = open(conf.conf['onePerLine'], 'r')

#make header for R
print 'density'
for line in resultsFile:	
	print line.strip().split('\t')[12]
	
