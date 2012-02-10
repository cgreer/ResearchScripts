import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy
import random
import cgWig

def plotAllDegOverlap(inFile, chrom, strand, wigDir1, wigDir2, outDir, withIntrons = False, flipStrand = True):
    '''hela must be 2nd wigDir2 cuz strand flip'''

    oppStrand = strand
    if flipStrand:
        oppStrand = bioLibCG.switchStrand(strand)
   
    print 'loading Wigs', chrom, strand
    coord_value1 = cgWig.loadSingleWig(wigDir1, chrom, strand, 'ALL')
    coord_value2 = cgWig.loadSingleWig(wigDir2, chrom, oppStrand, 'ALL')


    f = open(inFile, 'r')
    for line in f:
        gName, dChrom, dStrand, exonStarts, exonEnds = line.strip().split('\t')
        if dChrom != chrom or dStrand != strand:
            continue
        exonStarts = [int(x) for x in exonStarts.split(',')]    
        exonEnds = [int(x) for x in exonEnds.split(',')]    
        print 'Plotting', gName

        #create the span info for boxplots (JUST EXONS!!!)
        exons = zip(exonStarts, exonEnds)
        introns = [(x[0] + 1, x[1] - 1) for x in zip(exonEnds[:-1], exonStarts[1:])]
        iLengths = [x[0] - x[1] + 1 for x in zip(exonStarts[1:], exonEnds[:-1])] 
        all = exons[:]
       
        if withIntrons:
            all.extend(introns) 
            
        all.sort()
        tSpan = [('exon', x) if x in exons else ('intron', x) for x in all] 


        #gather expression data
        c_v = {}
        c_v2 = {}
        for type, (eStart, eEnd) in tSpan:
            for i in range(eStart, eEnd + 1):
                if i in coord_value1:
                    c_v[i] = coord_value1[i]
                if i in coord_value2:
                    c_v2[i] = coord_value2[i]

        #intron displacement for ONLY EXONS
        if not withIntrons:
            iCumulativeLengths = [sum(iLengths[:x]) for x in range(1,len(introns) + 1)]
            for i, (eStart, eEnd) in enumerate(exons):
                if i == 0: continue
                dAmount = iCumulativeLengths[i - 1] 
                for j in range(eStart, eEnd + 1):
                    if j in c_v:
                        c_v[j - dAmount] = c_v[j]
                        del c_v[j]
                    if j in c_v2:
                        c_v2[j - dAmount] = c_v2[j]
                        del c_v2[j]

        #get overall max
        overMax = max([max(x) for x in [c_v.values(), c_v2.values()]])
        a, b = set(c_v.keys()), set(c_v2.keys())
        overlap = a.intersection(b)
        colors_a = ['r' if x in overlap else 'k' for x in sorted(a)]
        colors_b = ['r' if x in overlap else 'k' for x in sorted(b)]

        plotGrassTrack(c_v, [9, 15], manualMax = overMax, flip = False, colors = colors_a)
        plotGrassTrack(c_v2, [-3, 3], manualMax = overMax,flip = True, colors = colors_b)
        xStart = plotGeneTrack(tSpan, 0)

        #labels and axes
        plt.figtext(.05, .5, gName)
        plt.figtext(.05, .62, '0 -')
        plt.figtext(.05, .89, '%s -' % overMax)
        plt.figtext(.05, 1 - .62, '0 -')
        plt.figtext(.05, 1 - .89, '%s -' % overMax)
        plt.ylim(-3,15)
        frame1 = plt.gca()
        frame1.axes.get_yaxis().set_visible(False)
        if dStrand == '1':
            plt.title('Degradome Comparison (5-->3)')
        else:
            plt.title('Degradome Comparison (3-->5)')
        imgName = outDir + '/' + gName + '.degOverlapPlot.png'
        plt.savefig(imgName, bbox_inches='tight', pad_inches=1)
        #plt.show()
        plt.close('all')
    
    f.close()


def plotGrassTrack(coord_value, yRange, manualMax = None, flip = False, colors = 'k'):
    
    # make adjusted dictionary
    # scale to manualMax and yRange

    #scale the length of the peak
    if not manualMax:
        manualMax = max(coord_value.values())
    maxScale = yRange[1] - yRange[0]        
    adjustedCoord_value = dict( (i, (float(j)/manualMax) * maxScale ) for i,j in coord_value.items())    

    #plot lines
    coords = adjustedCoord_value.keys()
    coords.sort()
    if flip:
        variedPoints = [yRange[1] - x for x in adjustedCoord_value.values()]
        plt.vlines(coords, variedPoints, yRange[1], colors = colors)
    else:
        variedPoints = [yRange[0] + x for x in adjustedCoord_value.values()]
        plt.vlines(coords, yRange[0], variedPoints, colors = colors)


def plotGeneTrack(tSpan, track):
    type_width = dict( (x[0], x[1]) for x in [('utr', 2), ('exon', 3), ('intron', 1)] )
    boxLengths = [pair[1] - pair[0] for type, pair in tSpan]
    totalLength = sum(boxLengths)
    firstCoord = tSpan[0][1][0]

    coveredLength = 0
    trackStart = track * 14 + (2* type_width['exon'])
    for i, spanInfo in enumerate(tSpan):
        currentLength = boxLengths[i]
        width = type_width[spanInfo[0]]
        x1, x2 = firstCoord + coveredLength, firstCoord + coveredLength + currentLength
        y1, y2 = trackStart - width, trackStart + width 

        plt.fill([x1,x2,x2,x1], [y1, y1, y2, y2], 'b', fill = False)

        coveredLength += currentLength

    return trackStart


if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

