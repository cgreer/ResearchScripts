import cgConfig as c
import bioLibCG as cg
import cgSort

mConf = c.cgConfig('Main.conf')

smallPath = mConf.conf['smallPath']

smallPath = '/home/chrisgre/smallLibs/WIGS/zebrafish'
#grab everything - NOT WIG MERGES...
smallLibs = cg.recurseDir(smallPath,  end = '.wig')
smallLibs.extend(cg.recurseDir(smallPath, end = '.wig'))

for lib in smallLibs:
	
	print 'sorting', lib
	cgSort.wigSort(lib)
