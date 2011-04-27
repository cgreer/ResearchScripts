import cgLuckyCharms
import bsddb3
#import tc

class Field:

	def __init__(self, dataType, dataDefault):
		self.dataType = dataType
		self.dataDefault = dataDefault

class DBBind:

	def __init__(self, dbFN, attributeField):
		
		#open DB
                self._db =  bsddb3.db.DB()
		#self._db.set_cachesize(0,524288000) #500 MB cache for each attribute...
		self._db.set_cachesize(1,0) #500 MB cache for each attribute...
		self._db.open(dbFN, None, bsddb3.db.DB_HASH, bsddb3.db.DB_CREATE)
		self.casteToDBFxn = cgLuckyCharms.getCasteFunction(attributeField.dataType, False )
		self.casteFromDBFxn = cgLuckyCharms.getCasteFunction(attributeField.dataType, True)
		self.defaultValue = attributeField.dataDefault

        def __iter__(self):
                '''got to be a better way, maybe write a custom iterator for the db...'''
                return iter(self._db.keys())

	def set(self, key, value):
		self._db.put(key, self.casteToDBFxn(value)) 

	def __setitem__(self, key, value):
		self._db[key] = self.casteToDBFxn(value)
	
	def get(self, key, default = "notPassed"):
                
                ret = self._db.get(key)
                if not ret:
                        if default == "notPassed":
                                return self.defaultValue
                        else:
                                return default
                else:
                        return self.casteFromDBFxn(ret)

	def __getitem__(self, key):
		return self.casteFromDBFxn(self._db.get(key))

	def commit(self):
		self._db.flush()

	def close(self):
		self._db.close()

class Nexus:

	def __init__(self, dataDir, dataClass):
		self.dataDir = dataDir
		self.dataClass = dataClass
	
	def bindAttribute(self, attributeName):
		
		#open database with create option hash
		setattr(self, attributeName, DBBind(self.dataDir + '/' + attributeName + '.db', getattr(self.dataClass, attributeName)))
		#bind it to the object
	
	def load(self, attNames):
		
		for attName in attNames:
			self.bindAttribute(attName)

	def close(self):
		pass	
