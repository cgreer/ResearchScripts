import bioLibCG
from cgNexusFlat import Field
                        

class OriginRNA:

        tcc = Field('string', None, 1)
        tccs = Field('stringList', list(), 2)
        sequence = Field('string', None, 3)
        eLevel = Field('int', 0, 4)
        context2 = Field('string', '.', 5)
        transcriptOverlap = Field('bool', False, 6)
        targets = Field('intList', list(), 7)
        entropy = Field('float', 0.0, 8)
        filteredTargets = Field('intList', list(), 9)
        endContigLength = Field('int', None, 10)
        totalContigLength = Field('int', None, 11)
        sequenceDuplicate = Field('bool', False, 12)
        targetDuplicate = Field('bool', False, 13)
        snr = Field('float', -1.0, 14)
        passedFilter = Field('bool', False, 15)
        avgNumSimulationTargets = Field('float', None, 16)
        transcriptIDs = Field('intList', list(), 17)
        transcriptContexts = Field('stringList', list(), 18)
        transcriptTypes = Field('stringList', list(), 19)
        transcriptType = Field('string', '.', 20)
        repeatStatus = Field('bool', False, 21)
        gScore = Field('int', None, 22)
        context = Field('string', '.', 23)
        numSignificantSequences = Field('int', 0, 24)
        avgNumSS = Field('float', -1.0, 25)
        snrSS = Field('float', 0.0, 26)
        phastScores = Field('floatList', list(), 27)
