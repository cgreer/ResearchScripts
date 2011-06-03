import bioLibCG
import cgNexusFlat

def checkMessy(tranFN):
        
        p = bioLibCG.cgPrint()                               
        f = open(tranFN, 'r')
        a = 0
        b = 0
        c = 0
        d = 0
        e = 0
        for line in f:
                ls = line.strip().split('\t')
                chrom, strand = ls[1], bioLibCG.switchStrandFormat(ls[2])
                tStart, tEnd = int(ls[3]), int(ls[4]) - 1
                cStart, cEnd = int(ls[5]), int(ls[6]) - 1
                exonStarts = [int(x) for x in ls[8][:-1].split(',')]
                exonEnds = [int(x) - 1 for x in ls[9][:-1].split(',')]
                exonPairs = zip(exonStarts, exonEnds)
                codingStatus = '_coding' in ls[15]
                tID = ls[0]

                #debug
                p.show = False

                intronPairs = []
                i = 0
                for pair in exonPairs:
                        if i == 0:
                                i += 1
                                continue
                        iStart = exonPairs[i -1][1] + 1
                        iEnd = exonPairs[i][0] - 1
                        intronPairs.append((iStart, iEnd))
                        i += 1

                #p.tell(tStart, tEnd, cStart, cEnd, exonPairs, intronPairs) 
                
                #take care of messy UTRs and assign utr ranges
                #5UTR
                if strand == '1':
                        if cStart == tStart or cStart == tEnd + 1:
                                p.tell('5 is none')
                                b += 1
                                if codingStatus:
                                        d += 1
                                range5 = ()
                        else:
                                range5 = (tStart, cStart - 1)
                else:
                        if cEnd + 1 == tStart or cEnd + 1 == tEnd + 1:
                                p.tell('5 is none')
                                b += 1
                                if codingStatus:
                                        d += 1
                                range5 = ()
                        else:
                                range5 = (cEnd + 1, tEnd)

                
                #3UTR
                if strand == '1':
                        if cEnd + 1 == tStart or cEnd + 1 == tEnd + 1:
                                p.tell('3 is none')
                                c += 1 
                                if codingStatus:
                                        e += 1
                                range3 = ()
                        else:
                                range3 = (cEnd + 1, tEnd)
                else:
                        if cStart == tStart or cStart == tEnd + 1:
                                p.tell('3 is none')
                                c += 1
                                if codingStatus:
                                        e += 1
                                range3 = ()
                        else:
                                range3 = (tStart, cStart - 1)

                a += 1

        print a, b, c, d, e                

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
