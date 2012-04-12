from cgNexusFlat import Field
import cgDegPeak

class cgTest:

	geneName = Field('string', '.', 1)                                                                             
	geneLength = Field('int', -1, 2)
	geneTargets = Field('intList', list(), 3)
