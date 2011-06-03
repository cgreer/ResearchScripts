import bioLibCG
import cgNexusFlat

def addOne(fN):
        
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                chrom, strand, start, end = bioLibCG.tccSplit(ls[0])
                start += 1
                end += 1

                print bioLibCG.makeTcc(chrom,strand,start,end)

def checkOverlaps(fn1, fn2):

        set1 = set()
        set2 = set()

        f = open(fn1, 'r')
        for line in f:
                ls = line.strip().split('\t')
                set1.add(ls[0])
        f.close()                
        print 'set 1: ', len(set1)

        f = open(fn2, 'r')
        for line in f:
                ls = line.strip().split('\t')
                set2.add(ls[0])
        f.close()                
        print 'set 2: ', len(set2)

        print '1 - 2'
        for i in set1 - set2:
                print i

        print '2 - 1'
        for i in set2 - set1:
                print i


if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
