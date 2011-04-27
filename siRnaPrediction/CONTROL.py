import siRnaPredict
import cgSort

'''
si.predictSiRNA('small.peak.data', 'siDegradome.conf')
si.updateSmallExpression('small.degradome', 'siPeaks.conf')
si.exonOverlap('small.degradome.small')
'''

'''
cgSort.sortFileCG('small.degradome.small.degOverlap', 1, rev = True)
cgSort.sortFileCG('small.degradome.small.degOverlap', 2, rev = True)
'''

'''
f = open('/home/chrisgre/smallLibs/siRNA/small/SRR029124.fastq.clipped.mapped', 'r')


for i in range(0,100): 
        r = f.readline()
        print r
        read = siRnaPredict.readInfo(r)
        print read.chrom, read.strand, read.start, read.end
        print read.scores, read.mismatches
        print read.returnMismatchPositions()

'''

directory = '/home/chrisgre/smallLibs/siRNA/small/sortedChroms'
siRnaPredict.markMismatches('small.degradome.small.degOverlap.clean', directory)
