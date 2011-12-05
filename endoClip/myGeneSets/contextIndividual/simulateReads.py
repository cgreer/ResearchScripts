import bioLibCG
import random
import cgNexusFlat
from progressbar import ProgressBar,SimpleProgress,Percentage,Bar,Timer

def simulateReads(fN, outFN):
        
        print 'creating context sets'
        c_s_set = {}
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                tcc = ls[2]
                chrom, strand, start, end = bioLibCG.tccSplit(tcc)
                for i in range(start, end + 1):
                        c_s_set.setdefault(chrom, {}).setdefault(strand, set()).add(i)
        f.close()

        print 'swapping sets for lists'
        for c in c_s_set:
                for s in c_s_set[c]:
                        c_s_set[c][s] = list(c_s_set[c][s])

        rChroms = c_s_set.keys()
        rStrands = ['1', '-1']
        print 'simulating reads'
        outF = open(outFN, 'w')
        pbar = ProgressBar(widgets=['  ', SimpleProgress(), ' ', Timer(), ' ', Bar()], maxval=8650000).start()
        for i in xrange(8650000):
                pbar.update(i)

                rChrom = random.choice(rChroms)
                rStrand = random.choice(rStrands)
                rCoord = random.choice(c_s_set[rChrom][rStrand])

                if rStrand == '1':
                        rStrand = '+'
                else:
                        rStrand = '-'

                pString = ['read_%s' % i, rStrand, rChrom, rCoord]
                outF.write('\t'.join([str(x) for x in pString]) + '\n')
                
        outF.close()                




if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
