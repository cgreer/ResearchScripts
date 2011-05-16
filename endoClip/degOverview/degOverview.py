import degPeak
import math
import cgWig
import cgNexusFlat
import bioLibCG
import gZipEntropy

def getEntropy(sequence):

        l = len(sequence)
        #get frequency of every nucleotide
        freqDict = {}
        for i in ['A', 'C', 'T', 'G']:
                freqDict[i] = sequence.count(i)/float(l)

        entropy = 0.0
        for i in ['A', 'C', 'T', 'G']:
                f = freqDict[i]
                if f == 0.0: continue
                entropy += (-1)*f*(math.log(f, 2))

        return entropy

def updateRepeatStatus(oFN, wigDir, chrom, strand, rn = None, tn = None):

        #load oRNAs
        oNX = cgNexusFlat.Nexus(oFN, degPeak.degPeak)
        oNX.load(['repeatStatus', 'tcc'])
        
        #load wig file for chrom, strand
        coord_value = cgWig.loadSingleWig(wigDir, chrom, strand, 'REPEAT')


        for oID in oNX.repeatStatus:

                oChrom, oStrand, start, end = bioLibCG.tccSplit(oNX.tcc[oID])
                if oChrom != chrom or oStrand != strand:
                        continue

                oNX.repeatStatus[oID] = False
                for i in range(start, end + 1):
                        if i in coord_value:
                                oNX.repeatStatus[oID] = True
                                break


        oNX.save()

def updateContext(oFN, wigDir, chrom, strand, rn = None, tn = None):

        oNX = cgNexusFlat.Nexus(oFN, degPeak.degPeak)
        oNX.load(['context', 'tcc'], [rn, tn])

        print 'loading wig'
        coord_contexts = cgWig.loadSingleWigContext(wigDir, chrom, strand, 'context') 
        print 'done loading'

        ds = bioLibCG.dominantSpotter(['C_EXON', 'C_3UTR', 'C_5UTR', 'NC_EXON', 'NC_3UTR', 'NC_5UTR', 'C_INTRON', 'NC_INTRON', 'INTER']) 


        for oID in oNX.tcc:

                oChrom, oStrand, start, end = bioLibCG.tccSplit(oNX.tcc[oID])
                
                #deg wigs is AS to actual clipping site
                if oStrand == '1':
                        oStrand = '-1'
                else:
                        oStrand = '1'
        
                if oChrom == chrom and oStrand == strand:

                        contexts = coord_contexts.get(start, 'INTER').split(',')
                        oNX.context[oID] = ds.spotItem(contexts)
        
        oNX.save()

def updateEntropy(oFN, rn = None, tn = None):
	
        oNX = cgNexusFlat.Nexus(oFN, degPeak.degPeak)
        oNX.load(['entropy', 'sequence'], [rn, tn])
        
        for oID in oNX.entropy:

                oNX.entropy[oID] = getEntropy(oNX.sequence[oID]) 

        oNX.save()

def updateSequence(oFN, seqFN, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, degPeak.degPeak)
        oNX.load(['sequence'], [rn, tn])

        f = open(seqFN, 'r')
        i = 0
        for line in f:
                ls = line.strip().split('\t')
                seq = ls[0]
                oNX.sequence[i] = seq
                i += 1

        oNX.save()

def updateTcc(oFN, peakFN, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, degPeak.degPeak)
        oNX.load(['tcc'], [rn, tn])

        f = open(peakFN, 'r')
        i = 0
        for line in f:
                ls = line.strip().split('\t')
                seq = ls[0]
                oNX.tcc[i] = seq
                i += 1

        oNX.save()

def updateELevel(oFN, wigDir, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, degPeak.degPeak)
        oNX.load(['tcc', 'eLevel'], [rn, tn])
        
        wigDict = cgWig.loadWigDict(wigDir)
        
        for oID in oNX.eLevel:
               
                
                expandTcc = bioLibCG.expandTcc(oNX.tcc[oID], 3)
                coord_value = cgWig.getExpressionProfile(expandTcc, wigDict)
                oNX.eLevel[oID] = max(coord_value.values())

        
        
        oNX.save()


def updateGScore(oFN, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, degPeak.degPeak)
        oNX.load(['gScore', 'sequence'], [rn, tn])

        for oID in oNX.gScore:

                oNX.gScore[oID] = gZipEntropy.gZipEntropy(oNX.sequence[oID])


        oNX.save()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
