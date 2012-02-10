import bioLibCG
import cgNexusFlat
import cgDL
from cgAutoCast import autocast
import GenomeFetch
from interval import IntervalSet, Interval
import cgWig
import matplotlib.pyplot as plt
import numpy as np
import math

class ASite:

    coord = cgNexusFlat.Field('string', '.', 1)                                                                             
    sequence = cgNexusFlat.Field('string', '.', 2)                                                                             

def updateSequence(dataFN, assembly = 'hg19', rn = None, tn = None):

    #load data
    NX = cgNexusFlat.Nexus(dataFN, ASite)
    NX.load(['coord', 'sequence'], [int(rn), int(tn)])

    #init GF
    gf = GenomeFetch.GenomeFetch(assembly)

    #update Sequence
    for id in NX.ids:
        chrom, strand, start, end = bioLibCG.tccSplit(NX.coord[id])
        NX.sequence[id] = gf.get_seq_from_to(chrom, start, end, strand)

    NX.save()

def locateASignals(dataFN, outFN, rn = None, tn = None):

    #load data
    NX = cgNexusFlat.Nexus(dataFN, ASite)
    NX.load(['coord', 'sequence'], [rn, tn])

    f = open(outFN, 'w')
    for id in NX.ids:
        chrom, strand, start, end = bioLibCG.tccSplit(NX.coord[id])
        if len(NX.sequence[id]) < 10: continue
        print NX.sequence[id], '\n'
        checkFrames = bioLibCG.returnFrames(NX.sequence[id], 6)
        for i, frame in enumerate(checkFrames):
            if frame == 'AATAAA':
                #assume 0-based...?
                siteStart, siteEnd = start + i, start + i + 5
                f.write('%s\n' % bioLibCG.makeTcc(chrom, strand, siteStart, siteEnd))
    f.close()

def plotBox(plotData):

    boxData = []
    f = open(plotData, 'r')
    for line in f:
        ls = line.strip().split('\t')
        boxData.append([int(x) for x in ls if x != '0'])
        print boxData[-1][:10]
    f.close()
    
    plt.boxplot(boxData, notch = 1, sym = '') 
    plt.show()

def getPlotData(aSites, wigDir, outFN):
    '''get box plot data from sites in degradome'''

    #load and init
    spreadRange = range(-200, 201) #200 +/- ... might want to check distance each AAUAAA is from each other
    relCoord_degVals = dict( (i, []) for i in spreadRange )
    
    for chrom in bioLibCG.humanChromosomes:
        for strand in ('1', '-1'):
            print chrom, strand
            coord_value = cgWig.loadSingleWig(wigDir, chrom, strand, 'ALL')
            f = open(aSites, 'r')
            for line in f:
                ls = line.strip().split('\t')
                ichrom, istrand, start, end = bioLibCG.tccSplit(ls[0])
                if ichrom != chrom or istrand != strand: continue

                for i in spreadRange:
                    degVal = coord_value.get(end + i, 0)
                    relCoord_degVals[i].append(degVal)
            f.close()
    
    #output box data
    #each row is a histogram of spread position (e.g., first row is -200)
    f = open(outFN, 'w')
    outLines = []
    for i in spreadRange:
        l = [str(x) for x in relCoord_degVals[i]]
        outLines.append('\t'.join(l) + '\n')
    f.writelines(outLines)
    f.close()

def plotSimplified(plotData):

    x, y, err = [], [], []
    f = open(plotData, 'r')
    for i, line in enumerate(f):
        x.append(i)
        Y, Err, a, b, c = line.strip().split('\t')
        Y, Err = float(Y), float(Err)
        y.append(Y)
        err.append(Err)
    f.close()

    plt.errorbar(x, y, err, fmt = 'o', ecolor='g')
    plt.show()
    

def simplifyPlotData(plotData, outFN):

    f = open(plotData, 'r')
    fOut = open(outFN, 'w')
    for line in f:
        data = [int(x) for x in line.strip().split('\t') if x != '0']
        data = np.array(data)
        n, avg, sd = data.size, data.mean(), data.std()
        error = sd / math.sqrt(n)
        fOut.write('%s\t%s\t%s\t%s\t%s\n' % (avg, error, n, sd, math.sqrt(n)))
    fOut.close()        
    f.close()
    
def collapseUTRs(utrFile, outFN):

    chrom_strand_iSet = {}
    f = open(utrFile, 'r')
    for line in f:
        ls = line.strip().split('\t')
        chrom, strand, start, end = bioLibCG.tccSplit(ls[0])
        chrom_strand_iSet.setdefault(chrom, {}).setdefault(strand, IntervalSet())
        chrom_strand_iSet[chrom][strand].add(Interval(start, end))
    f.close()

    fOut = open(outFN, 'w')
    for chrom in chrom_strand_iSet:
        for strand in chrom_strand_iSet[chrom]:
            iSet = chrom_strand_iSet[chrom][strand]
            for interv in iSet:
                try:
                    newTcc = bioLibCG.makeTcc(chrom, strand, interv.lower_bound, interv.upper_bound)
                    fOut.write('%s\n' % newTcc)
                except AttributeError:
                    print interv, 'was not added'

    fOut.close()

@autocast
def get3UTRFromTranscriptome(tranFN, outFN, wholeGene = False ):

        fOut = open(outFN, 'w')
        f = open(tranFN, 'r')
        for i, line in enumerate(f):
            ls = line.strip().split('\t')
            tChrom, tStrand = ls[1], bioLibCG.switchStrandFormat(ls[2])
            tStart, tEnd = int(ls[3]), int(ls[4]) - 1
            cStart, cEnd = int(ls[5]), int(ls[6]) - 1
            
            if wholeGene:
                utrTcc = bioLibCG.makeTcc(tChrom, tStrand, tStart, tEnd)
                fOut.write('%s\n' % utrTcc) 
                continue

            #5UTR
            if tStrand == '1':
                range5 = (tStart, cStart - 1)
            else:
                range5 = (cEnd + 1, tEnd)

            
            #3UTR
            if tStrand == '1':
                range3 = (cEnd + 1, tEnd)
            else:
                range3 = (tStart, cStart - 1)

            utrTcc = bioLibCG.makeTcc(tChrom, tStrand, range3[0], range3[1])
            fOut.write('%s\n' % utrTcc) 
        f.close()
        fOut.close()


if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

