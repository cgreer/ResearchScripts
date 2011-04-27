from cgNexusFlat import Field
import cgNexusFlat

class Seq:

	sequence = Field('string', '', 1)
	length = Field('int', 0, 2)
	filtered = Field('bool', False, 3)

def loadSeqs(seqFN):

	seqNX = cgNexusFlat.Nexus(seqFN, Seq)
	seqNX.load(['length', 'sequence'])

	print seqNX.length[100000], seqNX.sequence[100000]
	print 'done loading'


import sys
#loadSeqs(sys.argv[1])
