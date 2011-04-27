import fastQTypes, cgConfig
import bioLibCG as cg
import subprocess
import os, time
import clusterCheck

#init
mainConf = cgConfig.cgConfig('Main.conf')
metaFileNames = [mainConf.conf['metaFileName']]
wrapperShell = '/home/chrisgre/scripts/smallRNAProcessing/clipAdapter.sh'

def clipAdapter(fName, adapter = None, validate = False, oName = None, overwrite = True):
	
	#Check to see if the file exists:
	putativeN = fName.replace('.fastq','.clipped.fastq')
	if os.path.isfile(putativeN):
		if overwrite:
			print '  Overwriting file', putativeN
			os.remove(putativeN)
		else:
			print '  \nNOT OVERWRITING FILE', putativeN
			return 1
			 
	#If the adapter is none, try to find it in the small.meta file
	if adapter is None:
		baseFName = cg.getBaseFileName(fName) + '.counts'
		for metaFileName in metaFileNames:
			mFile = open(metaFileName, 'r')
			for line in mFile:
				fields = line.strip().split('\t')
				if baseFName == fields[0]:
					if fields[3] == 'NONE':
						print '  NO ADAPTER KNOWN FOR', fName
						return 1
					else:
						adapter = fields[3]
						print '  Using adapter', adapter, fName
			mFile.close()
	
	
	
	#Is it a valid fastq file?
	if validate: 
		pass
	
	#check the type of fastq file
	sangerType = False
	fType = fastQTypes.getFastQType(fName, quick = True)
	if fType == 'Sa':
		sangerType = True
	print '  Detected format:', fType, fName
	
	#Run it through clipper
	print 'Clipping file', fName
	if oName is None:
		oName = fName.replace('.fastq','.clipped.fastq')
	
	if sangerType:
		subprocess.Popen(['fastx_clipper', '-n', '-v', '-Q', '33', '-i', str(fName), '-a', str(adapter), '-o', str(oName)]).wait()
	else:
		subprocess.Popen(['fastx_clipper', '-n', '-v', '-i', str(fName), '-a', str(adapter), '-o', str(oName)]).wait()
	print '  DONE', fName
	
def clipAdapterInDirQ(dirName):
	'''The Q if for doing it on a cluster using qsub
	Every Q function has a corresponding shell script'''
	
	for file in cg.recurseDir(dirName, end = '.fastq'):
		
		#check if it isn't a clipped file:
		if 'clipped' in file:
			continue
			
		while True:
			#submit job if there are less than ten
			if clusterCheck.queryNumJobsQ('chrisgre') < 100:
				subprocess.Popen(['qsub', '-V', '-cwd', '-o', 'errors', '-e', 'errors', wrapperShell, file])
				time.sleep(.2) #give it time to update qstat
				break
			else:#wait 10 secs...
				time.sleep(10)

if __name__ == "__main__":
	import sys
	
	clipAdapter(str(sys.argv[1]))
