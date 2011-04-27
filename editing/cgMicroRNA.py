import bioLibCG
import cgSeqMod

class microRNA:

        def __init__(self):
                self.id = None
                self.name = None
                self.targetGenes = []
                self.seed = None
                self.sequence = None
                self.conservationStatus = None
                self.comSeed = None
                self.numBefore = 0
                self.numAfter = 0


def loadMicroRNAFromValidated(targetFN, sequenceFN):

        miNames_micros = {}

        f = open(targetFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                gName = ls[0]
                mNames = ls[1].split(',')
                
                for mName in mNames:
                        if mName not in miNames_micros:
                                m = microRNA()
                                m.name = mName
                                if gName not in m.targetGenes:
                                        m.targetGenes.append(gName)
                                miNames_micros[mName] = m
                        else:
                                m = miNames_micros[mName]
                                if gName not in m.targetGenes:
                                        m.targetGenes.append(gName)
        f.close()

        f = open(sequenceFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                miName = ls[0]
                miSequence = ls[1]

                if miName in miNames_micros:
                        m = miNames_micros[miName]
                        m.sequence = miSequence
                        m.seed = m.sequence[1:8]
                        m.comSeed = cgSeqMod.complementSequence(m.seed, True)

        f.close()

        print len(miNames_micros.values())
        count = 0
        for micro in miNames_micros.values():
                if micro.comSeed is None:
                        count += 1
        print count

        return miNames_micros
        
def loadMicroRNAFromTargetScan(fN, species):

        f = open(fN, 'r')
        f.readline() #header

        microRNAs = []
        for line in f:
                m = microRNA()
                ls = line.strip().split('\t')

                m.id = ls[3]
                #check species, there are also codes if you want to do it that way
                if not species in m.id:
                        continue
                
                m.seed = ls[1]
                m.sequence = ls[4]
                m.conservationStatus = int(ls[5])

                microRNAs.append(m)
        f.close()

        return microRNAs


