import bioLibCG
import cgNexusFlat
import cgDegPeak
import compareData

def countWithBins(dFN, binDir, type = 'INTRON'):

        print 'loading degradome'
        dNX = cgNexusFlat.Nexus(dFN, cgDegPeak.Peak)
        dNX.load(['tcc'])

        
        print 'loading bins'
        bins = []
        for i in range(50):
                bins.append([])

        for chrom in bioLibCG.humanChromosomes:
                for strand in ('1', '-1'):
                        f = open(binDir + '/%s.%s.%s.bins' % (type, chrom, strand), 'r')
                        for line in f:
                                ls = line.strip().split('\t')
                                tccs = ls[1:51]
                                for i in range(0,50):
                                        bins[i].append(tccs[i])

        #collect dTtcs in list
        dTccs = []
        for dID in dNX.tcc:
                tcc = dNX.tcc[dID]
                c, s, st, en = bioLibCG.tccSplit(tcc)
                if s == '1':
                        s = '-1'
                        en = st
                else:
                        s = '1'
                        st = en
                dTccs.append(bioLibCG.makeTcc(c,s,st,en))

        print len(dTccs), len(bins[0])
        
        for i in range(0, 50):
                print i
                overlaps = compareData.compareTwoTcc(dTccs, bins[i], 1)
                print len(overlaps)
                print overlaps


def countWithBinsSetReads(readFN, binDir, type = 'INTRON'):

        numBins = 100
       
        #create bin sets
        c_s_bin_set = {}
        for chrom in bioLibCG.humanChromosomes:
                for strand in ('1', '-1'):
                        #initialize data structure
                        for i in range(0, numBins):
                                c_s_bin_set.setdefault(chrom, {}).setdefault(strand, {})[i] = set()
                        f = open(binDir + '/%s.%s.%s.bins' % (type, chrom, strand), 'r')
                        for line in f:
                                ls = line.strip().split('\t')
                                tccs = ls[1:numBins + 1]
                                for i in range(0,numBins):
                                        ch, st, sta, end = bioLibCG.tccSplit(tccs[i])
                                        for j in range(sta, end + 1):
                                                c_s_bin_set[chrom][strand][i].add(j)

        print 'creating read sets'
        #creat read set
        c_s_set = {}
        for chrom in bioLibCG.humanChromosomes:
                c_s_set[chrom] = {}
                for strand in ('1', '-1'):
                        c_s_set[chrom][strand] = set()

        f = open(readFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                strand, chrom, start = ls[1:4]
                if chrom not in bioLibCG.humanChromosomes: continue
                start = int(start)
                if strand == '+':
                        strand = '-1'
                        start += 20
                else:
                        strand = '1'

                c_s_set[chrom][strand].add(start)
                

        print 'counting'
        #make bCounts
        binCounts = [0] * numBins

        #count for each bin
        for i in range(0, numBins):
                for chrom in bioLibCG.humanChromosomes:
                        for strand in ('1', '-1'):
                                for j in c_s_set[chrom][strand]:
                                        if j in c_s_bin_set[chrom][strand][i]:
                                                binCounts[i] += 1
                
                print '%s\t%s' % (i, binCounts[i])


def countWithBinsSet(dFN, binDir, type = 'INTRON'):

        dNX = cgNexusFlat.Nexus(dFN, cgDegPeak.Peak)
        dNX.load(['tcc'])

        
        numBins = 1
        
        c_s_bin_set = {}

        for chrom in bioLibCG.humanChromosomes:
                for strand in ('1', '-1'):
                        #initialize data structure
                        for i in range(0, numBins):
                                c_s_bin_set.setdefault(chrom, {}).setdefault(strand, {})[i] = set()
                        f = open(binDir + '/%s.%s.%s.bins' % (type, chrom, strand), 'r')
                        for line in f:
                                ls = line.strip().split('\t')
                                tccs = ls[1:numBins + 1]
                                for i in range(0,numBins):
                                        ch, st, sta, end = bioLibCG.tccSplit(tccs[i])
                                        for j in range(sta, end + 1):
                                                c_s_bin_set[chrom][strand][i].add(j)

        #collect dTtcs in list
        dTccs = []
        for dID in dNX.tcc:
                tcc = dNX.tcc[dID]
                c, s, st, en = bioLibCG.tccSplit(tcc)
                if s == '1':
                        s = '-1'
                        en = st
                else:
                        s = '1'
                        st = en
                dTccs.append(bioLibCG.makeTcc(c,s,st,en))

        #make bCounts
        binCounts = [0] * numBins

        #count for each bin
        for i in range(0, numBins):
                for dTcc in dTccs:
                        c, s, st, en = bioLibCG.tccSplit(dTcc)
                        for j in range(st, en + 1):
                                if j in c_s_bin_set[c][s][i]:
                                        binCounts[i] += 1
                print '%s\t%s' % (i, binCounts[i])


if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
