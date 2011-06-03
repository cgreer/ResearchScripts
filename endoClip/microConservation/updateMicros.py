import bioLibCG
import cgNexusFlat
import cgOriginRNAFlat

def updateTccAndSNR(mFN, alignmentFN):
       
        mNX = cgNexusFlat.Nexus(mFN, cgOriginRNAFlat.OriginRNA)
        mNX.load(['tcc', 'snrSS'])

        f = open(alignmentFN, 'r')
        i = 0
        for line in f:
                ls = line.strip().split('\t')
                chrom = ls[2]
                strand = bioLibCG.switchStrandFormat(ls[1])
                start = int(ls[3]) + 1 # 1BASE conversion
                end = start + len(ls[4])

                newTcc = bioLibCG.makeTcc(chrom, strand, start, end)
                newSNR = 10.0

                mNX.tcc[i] = newTcc
                mNX.snrSS[i] = newSNR

                i += 1

        mNX.save()                

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
