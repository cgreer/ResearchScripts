import bioLibCG


def getBracketList(a):
        
        justCommas = a[1:-1].replace(' [', '|').replace('],', '').replace(' ', '').replace('[', '').replace(']', '').split('|')
        returnList = []
        for item in justCommas:
                returnList.append(item.split(','))

        return returnList                

def convertPsuedo(fN, nameFN):
        
        ensID_gID = {}
        f = open(nameFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                if ls[1] == '':
                        continue
                 
                ensID_gID[ls[0]] = ls[1]

        f = open(fN, 'r')
        f.readline()
        for line in f:
                ls = line.strip().split('\t')
                '''
                for i, text in enumerate(ls):
                        print i, text
                '''
                tID = ls[0]
                chr = 'chr' + ls[1]
                strand = ls[4]
                tss, tse = ls[2], ls[3]
                css, cse = tss, tss
                exons = getBracketList(ls[18])
                numExons = len(exons)
                exonStarts = ','.join([x[0] for x in exons])
                exonEnds = ','.join([x[1] for x in exons])
                geneName = ls[8]
                geneName = ensID_gID.get(geneName, geneName)
                stat5 = 'none'
                stat3 = 'none'
                tCoding = 'pseudogene_noncoding'
                unused = 'none'
                gCoding = 'pseudogene_noncoding'


                print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (tID, chr, strand, tss, tse, css, cse, numExons, exonStarts, exonEnds, geneName, stat5, stat3, tCoding, unused, gCoding)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(convertPsuedo, sys.argv)
