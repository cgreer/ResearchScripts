from interval import IntervalSet, Interval
import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast

def spacerDistData(tranFN, outFN):
    '''chr strand tranStart tranEnd'''

    chrom_length = bioLibCG.returnChromLengthDict('hg19')

    chrom_strand_iSet = {}
    for chrom in chrom_length:
        for strand in ('+', '-'):
            chrom_strand_iSet.setdefault(chrom, {}).setdefault(strand, IntervalSet())

    print 'making intervals'
    f = open(tranFN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        tranStart, tranEnd = int(ls[3]), int(ls[4])
        strand = ls[2]
        chrom = ls[1]

        chrom_strand_iSet[chrom][strand].add(Interval(tranStart, tranEnd))

    f.close()
    
    spacerData = []
    print 'creating spacer data'
    for chrom in chrom_strand_iSet:
        for strand in chrom_strand_iSet[chrom]:
            iSet = chrom_strand_iSet[chrom][strand]
            for i, interv in enumerate(iSet):
                if interv == iSet[-1]: break
                nextInterv = iSet[i + 1]
                seperation = nextInterv.lower_bound - interv.upper_bound
                spacerData.append(seperation)

    f = open(outFN, 'w')
    outLines = [str(x) + '\n' for x in spacerData]
    f.writelines(outLines)
    f.close()
     
    

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

