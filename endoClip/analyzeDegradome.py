import bioLibCG
import cgOriginRNAFlat
import math
import cgNexusFlat
import cgDegPeak
import cgIContext
import matplotlib.pyplot as plt

def pieFractions(countList):
        
        s = sum(countList)
        fracs = [float(x)/s for x in countList]
        return fracs

def degContextPie(degFN, iContextDir):

        #load degNX
        dNX = cgNexusFlat.Nexus(degFN, cgDegPeak.Peak)
        dNX.load(['tcc', 'iContexts', 'eLevel'])

        #load all context NXs
        print 'loading NXs'
        chrom_strand_NX = {}
        for chrom in bioLibCG.humanChromosomes:
                for strand in ('1', '-1'):
                        
                        iFN = iContextDir + '/iContext.%s.%s' % (chrom, strand) 
                        NX = cgNexusFlat.Nexus(iFN, cgIContext.IContext)
                        NX.load(['type', 'tcc'])

                        chrom_strand_NX.setdefault(chrom, {})[strand] = NX
        print 'done...'
       



        context_count = {}
        for dID in dNX.tcc:

                chrom, strand, start, end = bioLibCG.tccSplit(dNX.tcc[dID])

                if strand == '1':
                        strand = '-1'
                elif strand == '-1':
                        strand = '1'
                else:
                        print 'WHAT?!'
                
                consForPeak = set()
                for iConID in dNX.iContexts[dID]:

                        if iConID == -1:
                                consForPeak.add('Intergenic')
                                continue

                        conType = chrom_strand_NX[chrom][strand].type[iConID]
                        if conType == 'NC_INTRON':
                                consForPeak.add('Noncoding Intron')
                        elif conType == 'C_INTRON':
                                consForPeak.add('Coding Intron')
                        elif conType == 'C_EXON':
                                consForPeak.add('Coding Exon')
                        elif conType == 'NC_EXON':
                                consForPeak.add('Noncoding Exon')
                        elif conType == 'C_3UTR':
                                consForPeak.add('3\' UTR')
                        elif conType == 'C_5UTR':
                                consForPeak.add('5\' UTR')
                        elif conType == 'NC_5UTR':
                                consForPeak.add('Noncoding Exon')
                        elif conType == 'NC_3UTR':
                                consForPeak.add('Noncoding Exon')
                        else:
                                print 'what are you?', conType

                lCons = len(consForPeak)
                for con in consForPeak:
                        context_count[con] = context_count.get(con, 0) + 1.00/lCons


                        

        labels = sorted(context_count.keys())
        counts = [int(context_count[x]) for x in labels]
        fracs = pieFractions(counts)

        #add numbers to labels
        labels = ['%s (%s)' % (x, int(context_count[x])) for x in labels]
        
        #plot
        plt.title('Degradome Peak Context')
        plt.pie(fracs, labels=labels, colors = ('Gold', 'GoldenRod', 'MediumBlue', 'CornFlowerBlue', 'DarkRed', 'Teal', 'DarkSlateGray'), shadow = True)
        plt.show()

def smallContextPie(oFN, title = 'Small RNA Peak Context'):

        #load degNX
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['context'])


        context_count = {}
        for oID in oNX.context:

                consForPeak = set()
                conType = oNX.context[oID]
                if conType == 'NC_INTRON':
                        consForPeak.add('Noncoding Intron')
                elif conType == 'C_INTRON':
                        consForPeak.add('Coding Intron')
                elif conType == 'INTER':
                        consForPeak.add('Intergenic')
                elif conType == 'C_EXON':
                        consForPeak.add('Coding Exon')
                elif conType == 'NC_EXON':
                        consForPeak.add('Noncoding Exon')
                elif conType == 'C_3UTR':
                        consForPeak.add('3\' UTR')
                elif conType == 'C_5UTR':
                        consForPeak.add('5\' UTR')
                elif conType == 'NC_5UTR':
                        consForPeak.add('Noncoding Exon')
                elif conType == 'NC_3UTR':
                        consForPeak.add('Noncoding Exon')
                else:
                        print 'what are you?', conType
                
                for con in consForPeak:
                        context_count[con] = context_count.get(con, 0) + 1


                        

        labels = sorted(context_count.keys())
        counts = [int(context_count[x]) for x in labels]
        fracs = pieFractions(counts)

        #add numbers to labels
        labels = ['%s (%s)' % (x, int(context_count[x])) for x in labels]
        
        #plot
        plt.title(title)
        plt.pie(fracs, labels=labels, colors = ('Gold', 'GoldenRod', 'MediumBlue', 'CornFlowerBlue', 'DarkRed', 'Teal', 'DarkSlateGray'), shadow = True)
        plt.show()

def eLevelHistogramSmall(degFN):
        
        #load degNX
        dNX = cgNexusFlat.Nexus(degFN, cgOriginRNAFlat.OriginRNA)
        dNX.load(['eLevel'])

     
        histValues = []
        for oID in dNX.eLevel:
                histValues.append(dNX.eLevel[oID])

        histVals = [math.log(x, 10) for x in histValues]
        plt.hist(histVals, 50, facecolor = 'Gold')
        type = 'Small Peaks'
        plt.title('Expression Level for %s' % type)
        plt.xlabel('log(Expression Level)')
        plt.ylabel('Number of %s' % type)


        plt.show()

def eLevelHistogram(degFN):
        
        #load degNX
        dNX = cgNexusFlat.Nexus(degFN, cgDegPeak.Peak)
        dNX.load(['tcc', 'iContexts', 'eLevel'])

     
        histValues = []
        for oID in dNX.eLevel:
                histValues.append(dNX.eLevel[oID])

        histVals = [math.log(x, 10) for x in histValues]
        plt.hist(histVals, 50, facecolor = 'Gold')
        type = 'Degradome Peaks'
        plt.title('Expression Level for %s' % type)
        plt.xlabel('log(Expression Level)')
        plt.ylabel('Number of %s' % type)


        plt.show()

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
