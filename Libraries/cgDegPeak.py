import bioLibCG
from cgNexusFlat import Field
                        

class Peak:

        tcc = Field('string', '.', 1)
        sequence = Field('string', '.', 2)
        eLevel = Field('int', 0, 3)
        tOverlap = Field('bool', False, 4)
        context = Field('string', '.', 5)
        repeatStatus = Field('bool', False, 6)
        gSequence = Field('string', '.', 7)
        gScore = Field('int', 0, 8)
        iContexts = Field('intList', list(), 11)
        repeatCount = Field('int', 0, 12)
        totalContig = Field('int', 0, 13)
