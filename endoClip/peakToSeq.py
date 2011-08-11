import bioLibCG
import GenomeFetch

def peakToSeq(peakFN, extend, outFN, assembly):
        #extend is +25 for degradome and -6/-4 for oRNA
        extend = int(extend)
        gf = GenomeFetch.GenomeFetch(assembly)

        outF = open(outFN, 'w')
        f = open(peakFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                chrom, strand, start, end = bioLibCG.tccSplit(ls[0])
                start, end = start - extend, end + extend
                newTcc = bioLibCG.makeTcc(chrom, strand, start, end)
                outF.write(gf.getSequence(newTcc) + '\n')

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(peakToSeq, sys.argv)
