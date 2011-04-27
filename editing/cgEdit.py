import bioLibCG

class EditingSite:

        def __init__(self):
                self.ID = None
                self.chromosome = None
                self.coordinate = None
                self.tcc = None
                self.gene = None
                self.transcripts = []
                self.strand = None
                self.context = []
                self.microTargets = []
                self.eMicroTargets = []
                self.before = []
                self.after = []
                self.eRatio = None

def updateContextEditingSites(eList, contextFN):
        '''Note: This overrules the gene/transcript from the load function'''

        #make eID_editsite
        eID_eSite = {}
        for eSite in eList:
                eID_eSite[eSite.ID] = eSite

        f = open(contextFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                eID = int(ls[0])
                transcript = ls[2]
                gene = ls[1]
                tType = ls[3]

                eSite = eID_eSite[eID]
                if tType not in eSite.context: eSite.context.append(tType)
                if transcript not in eSite.transcripts: eSite.transcripts.append(transcript)
                eSite.gene = gene
                



def loadEditingSites(fN, nt = 'A'):
        '''Using our labs format, load the editing site into a list'''

        cBases = {'A':'T', 'T':'A', 'G':'C', 'C':'G'} 

        f = open(fN, 'r')
        f.readline() #header

        eList = []
        for line in f:
                ls = line.strip().split('\t')

                e = EditingSite()
                e.chromosome = ls[0]
                e.coordinate = int(ls[1])
                e.gene = ls[3]
                e.eRatio = ls[6]
                refBase = ls[2]
                cBase = cBases[nt]        
                if refBase == cBase:
                        e.strand = '-1'
                else:
                        e.strand = '1'
                
                e.tcc = bioLibCG.makeTcc(e.chromosome, e.strand, e.coordinate, e.coordinate)

                e.ID = int(ls[13])

                eList.append(e)
        f.close()

        return eList
