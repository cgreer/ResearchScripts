import fastQTypes
import subprocess, time, os
import cgConfig
import bioLibCG as cg
import clusterCheck

#init
mainConf = cgConfig.cgConfig('Main.conf')
metaFileNames = [mainConf.conf['metaFileName']]

def mapFastQ(fName, organism):
	
	indexFileHuman = '/home/chrisgre/indexes/bowtie/hg19'
	indexFileMouse = '/home/chrisgre/indexes/bowtie/mm9'
	indexFileZebrafish = '/home/chrisgre/indexes/bowtie/danRer6'
	
	if organism == 'human':
		indexName = indexFileHuman
	elif organism == 'mouse':
		indexName = indexFileMouse
	elif organism == 'zebrafish':
		indexName = indexFileZebrafish
		
	
	outName = fName + '.mapped'
	
		
	logFile = open(mainConf.conf['outLog'] + cg.getBaseFileName(fName), 'w')
	errorFile = open(mainConf.conf['errorLog'] + cg.getBaseFileName(fName), 'w')
	
	if fastQTypes.getFastQType(fName, quick = True) == 'Sa':
		print 'Mapping with 33 phred offset'
		subprocess.Popen(['bowtie', '--phred33-quals', '-k', '20', '-m', '20', '-p', '1', indexName, fName, outName], stdout=logFile, stderr=errorFile).wait()
	else:
		print 'Mapping with 64 phred offset'
		subprocess.Popen(['bowtie', '--phred64-quals', '-k', '20', '-m', '20', '-p', '1', indexName, fName, outName], stdout=logFile, stderr=errorFile).wait()
	
	logFile.close()
	errorFile.close()
	
def mapFastQInDirQ(dirName, overwrite = True):
	'''Every Q function has a corresponding shell script'''
	wrapperShell = '/home/chrisgre/scripts/mapping/mapFastQ.sh'
	
	
	for file in cg.recurseDir(dirName, end = 'clipped.fastq'):
		print file
		
		putativeN = file.replace('.clipped.fastq','.clipped.fastq.mapped')
		if os.path.isfile(putativeN):
			if overwrite:
				print '  Overwriting file', putativeN
				os.remove(putativeN)
			else:
				print '  \nNOT OVERWRITING FILE', putativeN
				continue
				
		#check if mouse or human
		baseFName = cg.getBaseFileName(file, naked = True)
		org = 'None'
		for metaFileName in metaFileNames:
			mFile = open(metaFileName, 'r')
			for line in mFile:
				fields = line.strip().split('\t')
				if baseFName == fields[0]:
					if fields[2] == 'NONE':
						print '  NO ORG KNOWN FOR', file
						continue
					else:
						org = fields[2]
						print '  USING ORG', org, file
			mFile.close()
			
		#check if there is an organism, must check due to files not in metaFile
		if org == 'None':
			print '  NO org', file
			continue
			
		while True:
			#submit job if there are less than ten
			if clusterCheck.queryNumJobsQ('chrisgre') < 40:
				#subprocess.Popen(['qsub', '-q', 'xiao', '-l', 'mem=4G', '-V', '-o', mainConf.conf['outLog'], '-e', mainConf.conf['errorLog'], wrapperShell, file, org ])
				subprocess.Popen(['qsub', '-l', 'mem=4G', '-V', '-o', mainConf.conf['outLog'], '-e', mainConf.conf['errorLog'], wrapperShell, file, org ])
				time.sleep(.2) #give it time to update qstat
				break
			else:#wait 10 secs...
				time.sleep(20)

if __name__ == "__main__":
	import sys
	
	mapFastQ(sys.argv[1], sys.argv[2])
