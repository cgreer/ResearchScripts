import bioLibCG
import cgNexusFlat
import cgOriginRNAFlat
import cgAlignmentFlat
import cgDegPeak
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats
import math
import numpy
from cgAutoCast import autocast


@autocast
def quickECDF(fN, columns, properties = '', hist = False, numBins = 100):
    '''hist only works when there is one column...
    multi columns works only for the same file...'''

    try:
        columns = [int(x) for x in columns.split(',')]
        properties = [x for x in properties.split(',')]
    except:
        columns = [int(columns)]
        properties = [properties]

    #collect all values
    column_values = {} 
    f = open(fN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        for column in columns:
            try:
                column_values.setdefault(column, []).append(float(ls[column]))
            except IndexError:
                pass
        
    sortedKeys = column_values.keys()
    sortedKeys.sort()
    
    if hist:
        plt.subplot(len(columns), 1, 1)

    #plot all ecdfs/hist
    for i, col in enumerate(sortedKeys):

        if not hist:
            plt.hist(column_values[col], bins = numBins, cumulative = True,
                    histtype = 'step', normed = True, label = properties[i])
        else:
            plt.subplot(len(columns), 1, i)
            plt.hist(column_values[col], alpha = .3, bins = numBins, label = properties[i])
            plt.xlabel('hist for %s' % properties[i])
            plt.legend()


    if not hist:
        plt.xlabel('ecdf for %s' % properties)
        plt.legend()
    plt.show()

@autocast
def quickBoxPlot(fN, columns, properties = '', logTransform = False):
    '''hist only works when there is one column...
    multi columns works only for the same file...'''

    try:
        columns = [int(x) for x in columns.split(',')]
        properties = [x for x in properties.split(',')]
    except:
        columns = [int(columns)]
        properties = [properties]

    #collect all values
    column_values = {} 
    f = open(fN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        for column in columns:
            if logTransform:
                theVal = float(ls[column])
                if theVal <= 0:
                    print 'skipping log value <= 0'
                    continue
                theVal = math.log(float(theVal), 10)
            else:
                theVal = float(ls[column])
            column_values.setdefault(column, []).append(theVal)
        
    sortedKeys = column_values.keys()
    sortedKeys.sort()
    

    #plot all ecdfs/hist
    boxData = []
    for i, col in enumerate(sortedKeys):
        boxData.append(column_values[col])

    plt.boxplot(boxData, notch = 1) 
    plt.xlabel('ecdf for %s' % properties)
    plt.legend()
    plt.show()

@autocast
def quickSummary(fN, columns, properties = ''):
    '''hist only works when there is one column...
    multi columns works only for the same file...'''

    try:
        columns = [int(x) for x in columns.split(',')]
        properties = [x for x in properties.split(',')]
    except:
        columns = [int(columns)]
        properties = [properties]

    #collect all values
    column_values = {} 
    f = open(fN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        for column in columns:
            theVal = float(ls[column])
            column_values.setdefault(column, []).append(theVal)
        
    sortedKeys = column_values.keys()
    sortedKeys.sort()
    

    #plot all ecdfs/hist
    boxData = []
    for i, col in enumerate(sortedKeys):
        numPoints, (theMin, theMax), theMean, theVar, theSkew, theKurt = stats.describe(column_values[col])
        theMedian = numpy.median(column_values[col])
        theSD = numpy.std(column_values[col])
        cumFreqs = list(plt.hist(column_values[col], cumulative=True, bins = 10)[1])
        print 'Properties for', properties[i], ':'
        print ' # dPoints', numPoints
        print ' Min', theMin
        print ' Max', theMax
        print ' Mean', theMean
        print ' Median', theMedian
        print ' SD', theSD
        print ' Variance', theVar
        print ' Skew', theSkew
        print ' Kurt', theKurt
        print ' Percentiles', cumFreqs



@autocast
def quickPlot(fN, columns, properties = ''):
    '''hist only works when there is one column...
    multi columns works only for the same file...'''

    try:
        columns = [int(x) for x in columns.split(',')]
        properties = [x for x in properties.split(',')]
    except:
        columns = [int(columns)]
        properties = [properties]

    #collect all values
    column_values = {} 
    f = open(fN, 'r')
    for line in f:
        ls = line.strip().split()
        for column in columns:
            column_values.setdefault(column, []).append(float(ls[column]))
        
    sortedKeys = column_values.keys()
    sortedKeys.sort()
    
    #plot all ecdfs/hist
    for i, col in enumerate(sortedKeys):
        
        plt.plot(column_values[col], label = properties[i])

    plt.xlabel('plot for %s' % properties)
    plt.legend()
    plt.show()

@autocast
def quickPie(fN, column, property = ''):
    '''hist only works when there is one column...
    multi columns works only for the same file...'''

    type_count = {}
    f = open(fN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        theType = ls[column]
        type_count[theType] = type_count.get(theType, 0) + 1
    
    labels = []
    counts = []
    for theType, count in type_count.items():
        labels.append(theType)
        counts.append(count)

    #plot all ecdfs/hist
    plt.pie(counts, labels = [str(x) for x in labels])
    plt.title(property)
    plt.show()

@autocast
def quickPie(fN, column, property = ''):
    '''hist only works when there is one column...
    multi columns works only for the same file...'''

    type_count = {}
    f = open(fN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        theType = ls[column]
        type_count[theType] = type_count.get(theType, 0) + 1
    
    labels = []
    counts = []
    for theType, count in type_count.items():
        labels.append(theType)
        counts.append(count)

    #plot all ecdfs/hist
    plt.pie(counts, labels = [str(x) for x in labels])
    plt.title(property)
    plt.show()
        
def quickECDFContext(fN, column, property = ''):

        contexts = ['C_EXON', 'NC_EXON', 'C_INTRON', 'NC_INTRON', 'C_5UTR', 'C_3UTR']

        context_values = {}

        f = open(fN, 'r')
        for line in f:
            ls = line.strip().split('\t')

            context = ls[5]
            context_values.setdefault(context, []).append(int(ls[int(column)]))
            
        for context, values in context_values.items():
            plt.hist(values, bins = 1000000, cumulative = True, histtype = 'step', normed = True, label = context)


        plt.xlabel('ecdf for %s' % property)
        plt.legend()
        plt.show()

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])        
