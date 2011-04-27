import bioLibCG
import cgDB

                        

class OriginRNA:

        id = cgDB.Field('int', None)
        tcc = cgDB.Field('string', None)
        tccs = cgDB.Field('stringList', [])
        sequence = cgDB.Field('string', None)
        eLevel = cgDB.Field('int', 0)
        microOverlap = cgDB.Field('bool', False)
        transcriptOverlap = cgDB.Field('bool', False)
        targets = cgDB.Field('intList', [])
        entropy = cgDB.Field('float', 0.0)
        filteredTargets = cgDB.Field('intList', [])
        endContigLength = cgDB.Field('int', None)
        totalContigLength = cgDB.Field('int', None)
        sequenceDuplicate = cgDB.Field('bool', False)
        targetDuplicate = cgDB.Field('bool', False)
        snr = cgDB.Field('float', None)
        passedFilter = cgDB.Field('bool', False)
        avgNumSimulationTargets = cgDB.Field('float', None)
        transcriptIDs = cgDB.Field('intList', [])
        transcriptContexts = cgDB.Field('stringList', [])
        transcriptTypes = cgDB.Field('stringList', [])
        transcriptCodingTypes = cgDB.Field('stringList', [])

        def __init__(self, id):
                self.id = id 




def prettyPrint(oRNA, message = ''):
        print message, oRNA.id, oRNA.entropy, oRNA.filteredTargets, oRNA.endContigLength, oRNA.totalContigLength, oRNA.sequenceDuplicate, oRNA.avgNumSimulationTargets, oRNA.snr
