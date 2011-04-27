import cgNexusFlat
import cgSeqFlat
import bioLibCG



def updateLength(dataFN, runNumber = None, totalRuns = None):

        #c = {'sequence' : lambda x: x in set(['Burger'])}        
        c = {}
	seqNX = cgNexusFlat.Nexus(dataFN, cgSeqFlat.Seq)
        seqNX.load(['length', 'sequence'], [runNumber, totalRuns], conditions = c )

	for id, val in seqNX.length.iteritems():
		seqNX.length[id] = val + 1

        #print seqNX.sequence, seqNX.length
	seqNX.save()

if __name__ == '__main__':
	import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

