#This handles all configuration file stuff
configDir = '/home/chrisgre/configFiles/'


#get the current configuration:
currentFile = open(configDir + 'current.meta')
currentConfig = currentFile.readline().strip()
currentFile.close()

#currentConfig = 'HY.conf'
#currentConfig = 'test.conf'


class cgConfig:

	conf = {}
	configName = ''
	
	def __init__(self, configFilename = currentConfig):
		self.configName = configFilename
		self.returnConfigDict(self.configName)
		

	def returnConfigDict(self, config = configName):
		configFileName = configDir + config
		configFile = open(configFileName, 'r')
		
		for line in configFile:
			if line.startswith('#') or line.startswith('\'') or (line == '\n'):
				pass
			else:#parse config
				
				self.conf[line.strip().split('::')[0]] = line.strip().split('::')[1]
				
	def printConfig(self):
		print 'CONFIGURATION(%s)' % self.configName
		for key in self.conf:
			print key, self.conf[key]

def getConfig(cName):
	if not cName:
		return cgConfig()
	else:
		return cgConfig(cName)
		
