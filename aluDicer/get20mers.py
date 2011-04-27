import GenomeFetch
import bioLibCG

def get20mers(aluFN):
        
        gf = GenomeFetch.GenomeFetch('hg19')
        seq_count = {}
        f = open(aluFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                coord = ls[0]
                chrom, start, end = coord.split(':')[0], coord.split(':')[1].split('-')[0], coord.split(':')[1].split('-')[1]
                strand = bioLibCG.switchStrandFormat(ls[2])
                tcc = bioLibCG.makeTcc(chrom, strand, start, end)
                seq = gf.getSequence(tcc)
                frames = bioLibCG.returnFrames(seq, 20)
                if frames == 1:
                        continue
                for smallSeq in frames:
                        count = seq_count.get(smallSeq, 0)
                        seq_count[smallSeq] = count + 1

        
        for seq, count in seq_count.items():
                print '%s\t%s' % (seq, count)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(get20mers, sys.argv)
