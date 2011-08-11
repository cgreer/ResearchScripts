import bioLibCG
import cgNexusFlat
import GenomeFetch

        
def makeTranscriptome(tranFN, outFN):

        p = bioLibCG.cgPrint()                               
        p.show = False
        gf = GenomeFetch.GenomeFetch('hg19')
        

        fOut = open(outFN, 'w')
        f = open(tranFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                tChrom, tStrand = ls[1], bioLibCG.switchStrandFormat(ls[2])
                exonStarts = [int(x) + 1 for x in ls[8][:-1].split(',')]
                exonEnds = [int(x) for x in ls[9][:-1].split(',')]
                exonPairs = zip(exonStarts, exonEnds)
                tID = ls[0]
                gID = ls[10]

                seqList = []
                for eStart, eEnd in exonPairs:
                        tcc = bioLibCG.makeTcc(tChrom, tStrand, eStart, eEnd)
                        seqList.append(gf.getSequence(tcc))

                mRNA = ''.join(seqList)

                #reverse direction if negative strand
                if tStrand == '-1':
                        mRNA = mRNA[::-1]

                fOut.write('> %s:%s:%s\n' % (tID, gID, len(mRNA)))
                fOut.write(mRNA + '\n\n')
                        
        fOut.close()
        f.close()



if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
