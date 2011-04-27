import bioLibCG
import GenomeFetch

def getPeakSequences(peakFN, extend = 0):
        
        extend = int(extend)

        f = open(peakFN, 'r')
        peaks = [x.strip() for x in f]
        f.close()
        
        gf = GenomeFetch.GenomeFetch('hg19')


        for peak in peaks:
                chrom, strand, start, end = bioLibCG.tccSplit(peak)
                start = start - extend
                end = end + extend
                print gf.get_seq_from_to(chrom, start, end, strand)


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(getPeakSequences, sys.argv)

