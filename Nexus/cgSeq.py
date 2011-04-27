from cgNexus import Field
import cgNexus

class Seq:

	sequence = Field('string', '')
	length = Field('int', 0)


def loadSeqs(seqFN):

	seqNX = cgNexus.Nexus('seqTest', Seq)
	seqNX.load(['length', 'sequence'])

	f = open(seqFN, 'r')
	for line in f:
		ls = line.strip().split('\t')
		id, seq, length = ls[0], ls[1], int(ls[2])
		seqNX.length[id] = length
		seqNX.sequence[id] = seq

	seqNX.length.close()
	seqNX.sequence.close()

import sys
#loadSeqs(sys.argv[1])
