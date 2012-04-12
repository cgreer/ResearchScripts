import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast

def collectIDs2(fN, fN2, fN3):
    '''Used for getting the # repeat reads on target results'''

    idSet = set()
    f = open(fN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        
        for id in ls[9].split(','):
            idSet.add(id)

    idSetDeg = set()
    f = open(fN2, 'r')
    for line in f:
        ls = line.strip().split('\t')
        
        if ls[18] == "F": continue

        if ls[0] in idSet:
            idSetDeg.add(ls[11])

    chrom_strand_coord = getHitMap(list(idSetDeg))

    #print 'target Tccs', len(idSetDeg)
    #print chrom_strand_coord['chr1']['1']

    eSet = set()
    eDict = {}
    readNames = set()
    f = open(fN3, 'r')
    for line in f:
        ls = line.strip().split('\t')
        chrom, strand, start = ls[2], ls[1], int(ls[3])
        strand = bioLibCG.switchStrandFormat(strand)
        for i in range(start - 3, start + 3):
            try:
                if i in chrom_strand_coord[chrom][strand]:
                    readNames.add(ls[0])
                    break
            except KeyError:
                continue
    f.close()

    #now go back through and count the times the read appears
    readName_count = {}
    f = open(fN3, 'r')
    for line in f:
        ls = line.strip().split('\t')
        if ls[0] in readNames:
            readName_count[ls[0]] = readName_count.get(ls[0], 0) + 1

    for read, count in readName_count.iteritems():
        print '%s\t%s' % (read, count)



def getHitMap(tccList):
	
        #create hitmap of chrom and strand
	chrom_strand_coord = {} #format = chr: { strand : { coord : value 
	for tcc in tccList:
                lChrom, lStrand, start, end = bioLibCG.tccSplit(tcc)
		lStrand = str(lStrand)
		start = int(start)
		end = int(end)

                if lChrom in bioLibCG.acceptableChroms:
               
                    #wig for regular
                    for i in range(start - 5, end + 5 ):
                        chrom_strand_coord.setdefault(lChrom, {}).setdefault(lStrand, set()).add(i)
        
        return chrom_strand_coord
    

def collectIDs(fN, fN2, fN3):
    '''collect ids from column of first file 
    filter lines that have those ids from second file'''

    idSet = set()
    f = open(fN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        
        for id in ls[9].split(','):
            idSet.add(id)

    idSetDeg = set()
    f = open(fN2, 'r')
    for line in f:
        ls = line.strip().split('\t')
        
        if ls[0] in idSet:
            idSetDeg.add(ls[2])
    
    f = open(fN3, 'r')
    for line in f:
        ls = line.strip().split('\t')
        
        if ls[0] in idSetDeg:
            print line.strip()

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
