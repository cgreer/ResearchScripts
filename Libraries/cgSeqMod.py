import bioLibCG

def loadCodonMap(assembly):
        '''Eventually Ill implement a map loading fxn, but for now this will do'''

        map = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
              "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
              "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*",
              "UGU":"C", "UGC":"C", "UGA":"*", "UGG":"W",
              "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
              "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
              "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
              "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
              "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
              "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
              "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
              "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
              "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
              "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
              "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
              "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
        
        return map

def getCodonListFromRNA(seq, readingFrame = 0):
        '''Given mRNA sequence, get list of each codon'''
        codonList = [seq[3*i:3*i + 3] for i in range(len(seq)//3)]
        return codonList

def translateRNA(seq, codonMap):
        
        codonList = getCodonListFromRNA(seq)
        translation = ''.join([codonMap[x] for x in codonList])


        return translation
        
def complementSequence(seq, rnaFlag):
        rnaFlag = bool(rnaFlag)

        rcRNA = {'A':'U', 'U':'A', 'C':'G', 'G':'C'}
        rcDNA = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

        if rnaFlag:
                rcDict = rcRNA
        else:
                rcDict = rcDNA

        newSeq = []
        for c in seq:
                newSeq.append(rcDict[c])
        
        complement = ''.join(newSeq)
        
        return complement

def reverseComplementSequence(seq, rnaFlag):
        rnaFlag = bool(rnaFlag)

        rcRNA = {'A':'U', 'U':'A', 'C':'G', 'G':'C'}
        rcDNA = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

        if rnaFlag:
                rcDict = rcRNA
        else:
                rcDict = rcDNA

        newSeq = []
        for c in seq:
                newSeq.append(rcDict[c])
        
        newSeq.reverse()
        complement = ''.join(newSeq)
        
        return complement

def transcribeDNA():
        pass

if __name__ == "__main__":
        import sys
        #bioLibCG.submitArgs(fxn, sys.argv)
        seq = 'AUUGCGCGACGUGUACCUAGUAUGCG'
        map = loadCodonMap('hg19')
        print seq
        print getCodonListFromRNA(seq)
        print translateRNA(seq, map)

