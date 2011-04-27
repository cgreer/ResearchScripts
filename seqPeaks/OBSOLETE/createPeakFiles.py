#using the continuos blocks from in the small RNA lib file, identify all the peaks...
import bioLibCG as cg
import cgConfig as c

mConf = c.cgConfig('Main.conf')

fileNames = cg.recurseDir(mConf.conf['wigMouse'], start = 'cb.', end = '.tsv')

for fN in fileNames:
	file = open(fN, 'r')
	
	#get peak positions and average heights
	lines = [str(int(x.split('\t')[0]) + int(x.split('\t')[2])/2) + '\t' + str(x.strip().split('\t')[3]) for x in file]
	file.close()
	
	
	
	
	
	
	
	
	#output new file with new ending
	outFile = open(fN + '.peakdata', 'w')
	for line in lines:
		outFile.write('%s\n' % line)
