import subprocess
import cgConfig
#check how many jobs per user

def queryNumJobsQ(user):
	
	#init
	mainConf = cgConfig.cgConfig('Main.conf')
	scratchFileName = mainConf.conf['qstatScratch']
	
	#output qstat info to scratch file
	sFile = open(scratchFileName, 'w')
	subprocess.Popen(['qstat', '-u', user], stdout=sFile).wait()
	sFile.close()
	
	#parse scratchfile, check how many times usr name appears
	sFile = open(scratchFileName, 'r')
	userCount = 0
	for line in sFile:
		if user in line:
			userCount = userCount + 1
	sFile.close()
	
	return userCount


if __name__ == "__main__":
	import sys
	
	queryNumJobsQ(sys.argv[1])
	
	
