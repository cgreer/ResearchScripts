import bioLibCG
import cgSeqMod

map = cgSeqMod.loadCodonMap('hg19')

numN = 0
numS = 0

for codon in map:
        bAlter = cgSeqMod.translateRNA(codon, map)
        codonList = list(codon)
        for i,nt in enumerate(codonList):
                if nt == 'A':
                        codonList[i] = 'G'
                        newCodon = ''.join(codonList)
                        newAlter = cgSeqMod.translateRNA(newCodon, map)
                        if newAlter == bAlter:
                                print 'SYN', bAlter, newAlter, codon, newCodon
                                numS += 1
                        else:
                                print 'NON', bAlter, newAlter, codon, newCodon
                                numN += 1

                        codonList[i] = 'A'

print numN, numS


