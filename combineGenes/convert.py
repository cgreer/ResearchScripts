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

def convertLongNC(fN, outFN):
        
        fOut = open(outFN, 'w')
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                
                chrom = ls[1]
                tss = int(ls[2])
                tse = int(ls[3])

                tName = ls[4]
                strand = ls[6]

                bLengths = [int(x) for x in ls[11][:-1].split(',')]
                bStarts = [int(x) for x in ls[12][:-1].split(',')]

                cStart = tss
                cEnd = tss

                numBlocks = len(bLengths)

                exonStarts = [tss + x for x in bStarts]
                exonEnds = [x + y for (x,y) in zip(exonStarts, bLengths)]

                if strand == '-':
                        exonStarts, exonEnds = exonEnds, exonStarts

                pString = [tName,
                           chrom,
                           strand,
                           tss,
                           tse,
                           cStart,
                           cEnd,
                           str(numBlocks),
                           ','.join([str(x) for x in exonStarts]),
                           ','.join([str(x) for x in exonEnds]),
                           'None',
                           'none',
                           'none',
                           'longNC_noncoding',
                           'longNC',
                           'noncoding_gene']
                pString = '\t'.join([str(x) for x in pString]) + '\n'
                fOut.write(pString)

        f.close()
        fOut.close()


if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

