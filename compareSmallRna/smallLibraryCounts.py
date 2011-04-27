###This is the MAIN script to compare small RNA libraries to the mature microRNA sequences.  It first calls
##collectMatureFrames and then calls compareFinished to do the actual comparisons.

from bioLibCG import *
import os
from collectMF import *
from compareCounts import *
import cgConfig

#init
conf = cgConfig.cgConfig()
baseName = conf.conf['baseName']
merLength = conf.conf['merLength']
inDirectory = conf.conf['inDirectory']
smallPath = conf.conf['smallPath']

mainConf = cgConfig.cgConfig('Main.conf')
smallPath = mainConf.conf['smallPath']
merLength = int(merLength)

####Make directory for files to go into
filePath = 'out/%s' % baseName
if baseName in os.listdir('out/'):
	#Delete all contents
	for file in os.listdir(filePath):
		os.remove('%s/%s' % (filePath,file) )
else:
	os.mkdir(filePath)

###Get ids and Xmer frames
collectMatureFrames(baseName, merLength)

###Compare against small libraries
compareCounts(baseName, merLength, smallPath)
	
