import bioLibCG
import cgNexusFlat
import random
import math

def getTotalLength(fN, outFN):
       
        lowCount = 0
        mult = 0.00077921994118935773
        tLen = 0
        reads = []
        fOut = open(outFN, 'w')
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                chrom = ls[1]
                strand = ls[2]
                start = int(ls[3])
                end = int(ls[4])
                thisLen = int(end - start)
                tLen += thisLen
                numReads = thisLen * mult
                reads.append(numReads)
                
                for i in range(math.ceil(numReads)):
                        j = random.randint(start, end)              
                        pString = ['read', strand, chrom, j]
                        fOut.write('\t'.join([str(x) for x in pString]) + '\n')
                

        print tLen                
        print min(reads), max(reads), sum(reads)
        print lowCount

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
