import fastQTypes
import bioLibCG as cg
import subprocess

metaFileNames = ['/u/home8/gxxiao/apps/projects/small.rna.libs/small.meta']
wrapperShell = '/something'

def clipAdapter(fName, adapter = None, validate = False, oName = None):
	
	#If the adapter is none, try to find it in the small.meta file
	if adapter is None:
		baseFName = cg.getBaseFileName(fName) + '.counts'
		for metaFileName in metaFileNames:
			mFile = open(metaFileName, 'r')
			for line in mFile:
				fields = line.strip().split('\t')
				if baseFName == fields[0]:
					if fields[3] == 'NONE':
						print 'NO ADAPTER KNOWN FOR', fName
						return 1
					else:
						adapter = fields[3]
						print 'Using adapter', adapter
			mFile.close()
	
	
	
	#Is it a valid fastq file?
	if validate: 
		pass
	
	#check the type of fastq file
	sangerType = False
	fType = fastQTypes.getFastQType(fName, quick = True)
	if fType == 'Sa':
		sangerType = True
	print 'Detected format:', fType
	
	#Run it through clipper
	print 'Clipping file'
	if oName is None:
		oName = fName.replace('.fastq','.clipped.fastq')
	
	if sangerType:
		subprocess.Popen(['fastx_clipper', '-v', '-Q', '33', '-i', str(fName), '-a', str(adapter), '-o', str(oName)]).wait()
	else:
		subprocess.Popen(['fastx_clipper', '-v', '-i', str(fName), '-a', str(adapter), '-o', str(oName)]).wait()
	print 'DONE'
	
def clipAdapterInDirQ(dirName):
	'''The Q if for doing it on a cluster using qsub
	Every Q function has a corresponding shell script'''
	
	for file in cg.recurseDir(dirName, end = '.fastq'):
		subprocess.Popen([wrapperShell, file])
