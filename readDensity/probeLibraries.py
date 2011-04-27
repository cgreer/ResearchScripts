import fastQTypes
import subprocess, time, os
import cgConfig as c
import bioLibCG as cg
import clusterCheck
import stepVectorScan


def probe(tcc, conf = None):
	
	if not conf:
		mConf = c.cgConfig('Main.conf')
	smallPath = mConf.conf['smallPath']
	
	chrom, strand, start, end = cg.tccSplit(tcc)
	
	total = 0
	for lib in cg.recurseDir(smallPath, end = 'mapped.%s.wig' % strand):
		
		
		try:
			eLevels = stepVectorScan.scanVectorsFile(lib, [tcc])
		except:
			print lib, 'index failed'
			continue
			
		
		#find highest expression level
		highest = 0
		for coord in eLevels:
			if eLevels[coord] > highest:
				highest = eLevels[coord]
				
				
		if highest > 0:
			print lib, highest
			total += highest
			#print eLevels
		
	print total
		
if __name__ == "__main__":
	import sys
	if len(sys.argv) < 3:
		sys.argv.append(None)
		
	probe(sys.argv[1], sys.argv[2])
		



