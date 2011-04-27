#This library clusters a list of tccCoordinates
import compareData as compare
import bioLibCG as cg

def returnTccListRange(cluster):
	chrom = cluster[0].split(':')[0]
	strand = cluster[0].split(':')[1]
	starts = [int(x.split(':')[2]) for x in cluster]
	ends = [int(x.split(':')[3]) for x in cluster]
	
	starts.sort()
	ends.sort()
		
	return '%s:%s:%s:%s' % (chrom, strand, starts[0], ends[-1])

class cgClusters:
	clusters = []
	addedList = []
	
	def createClusters(self, connections, item, mode):
			if item in self.addedList:
				#print item, 'addedList', self.addedList
				return 0
			elif mode == "top":
				self.clusters.append([item])
				self.addedList.append(item) ##creates new cluster with the item already stored in it
				#print item, 'top', self.addedList, self.clusters
				for connectedItem in connections[item]:
					self.createClusters(connections, connectedItem, "neighbor")
			elif mode == "neighbor":
				self.clusters[-1].append(item) #add this item to the last cluster created
				self.addedList.append(item)
				#print item, 'neigh', self.addedList, self.clusters
				for connectedItem in connections[item]:
					self.createClusters(connections, connectedItem, "neighbor")
					
	def addTest(self):
		self.clusters.append('hello')
		
	def makeClusters(self, tccList):
		#init
		self.addedList = []
		self.clusters = []
		
		#Create Clusters
		connectionsDict = compare.makeConnectionsDict(tccList)
		
		for tcc in tccList:
			self.createClusters(connectionsDict, tcc, "top")
			
		for cluster in self.clusters:
			cluster = cg.sortTccList(cluster)
			
	def getClusters(self):
		return self.clusters
		

