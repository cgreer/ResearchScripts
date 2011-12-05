import bioLibCG
from cgNexusFlat import Field
                        

class miR:

        mID = Field('string', '.', 1)
        accID = Field('string', '.', 2)
        sequence = Field('string', '.', 3)
        clippedSeq = Field('string', '.', 4)
        polySeqs = Field('stringList', list(), 5)
