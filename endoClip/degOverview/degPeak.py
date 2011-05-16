import bioLibCG
from cgNexusFlat import Field
                        

class degPeak:

        tcc = Field('string', None, 1)
        sequence = Field('string', None, 2)
        eLevel = Field('int', 0, 3)
        entropy = Field('float', 0.0, 4)
        repeatStatus = Field('bool', False, 5)
        gScore = Field('int', None, 6)
        context = Field('string', '.', 7)
