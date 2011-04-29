import bioLibCG
from cgNexusFlat import Field
                        

class OriginRNA:

        tcc = Field('string', None, 1)
        tccs = Field('stringList', list(), 2)
        sequence = Field('string', None, 3)
        eLevel = Field('int', 0, 4)
        microOverlap = Field('bool', False, 5)
        transcriptOverlap = Field('bool', False, 6)
        targets = Field('intList', list(), 7)
        entropy = Field('float', 0.0, 8)
        filteredTargets = Field('intList', list(), 9)
        endContigLength = Field('int', None, 10)
        totalContigLength = Field('int', None, 11)
        sequenceDuplicate = Field('bool', False, 12)
        targetDuplicate = Field('bool', False, 13)
        snr = Field('float', None, 14)
        passedFilter = Field('bool', False, 15)
        avgNumSimulationTargets = Field('float', None, 16)
        transcriptIDs = Field('intList', list(), 17)
        transcriptContexts = Field('stringList', list(), 18)
        transcriptTypes = Field('stringList', list(), 19)
        transcriptCodingTypes = Field('stringList', list(), 20)
        repeatStatus = Field('bool', False, 21)

