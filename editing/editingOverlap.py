import bioLibCG
import compareData
import cgEdit

def overlapWithDegradome(dFN, eFN):

        eSites = cgEdit.loadEditingSites(eFN)

        degTccs = []
        f = open(dFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                chrom, strand, start, end = bioLibCG.tccSplit(ls[1])
                start = start - 3
                end = end + 3
                degTccs.append(bioLibCG.makeTcc(chrom,strand,start,end))
        print degTccs[0:5]
        eTccs = [eSite.tcc for eSite in eSites]
       
        overlaps = compareData.compareTwoTcc(eTccs, degTccs, 1)

        print len(overlaps)


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(overlapWithDegradome, sys.argv)
        
