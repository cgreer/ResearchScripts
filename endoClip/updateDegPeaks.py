import bioLibCG
import cgDegPeak
import cgWig
import cgNexusFlat
import gZipEntropy

def updateTcc(oFN, tccFN, rn = None, tn = None):

        oNX = cgNexusFlat.Nexus(oFN, cgDegPeak.Peak)
        oNX.load(['tcc'], [rn, tn])

        f = open(tccFN, 'r')
        i = 0
        for line in f:
                ls = line.strip().split('\t')
                tcc = ls[0]
                oNX.tcc[i] = tcc 
                i += 1

        oNX.save()
        
def updateSequence(oFN, seqFN, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, cgDegPeak.Peak)
        oNX.load(['sequence'], [rn, tn])

        f = open(seqFN, 'r')
        i = 0
        for line in f:
                ls = line.strip().split('\t')
                seq = ls[0]
                oNX.sequence[i] = seq
                i += 1

        oNX.save()

def updateELevel(oFN, wigDir, rn = None, tn = None):
        '''Dont need to do it by chromosome because it is small enough'''
        '''Also dont need to flip the strand because the wig is opposite as well'''

        oNX = cgNexusFlat.Nexus(oFN, cgDegPeak.Peak)
        oNX.load(['tcc', 'eLevel'], [rn, tn])
        
        wigDict = cgWig.loadWigDict(wigDir)
        
        for oID in oNX.eLevel:
                
                expandTcc = bioLibCG.expandTcc(oNX.tcc[oID], 3)
                coord_value = cgWig.getExpressionProfile(expandTcc, wigDict)
                oNX.eLevel[oID] = max(coord_value.values())

        
       
        oNX.save()

def updateTranscriptOverlap(oFN, wigDir, chrom, strand, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, cgDegPeak.Peak)
        oNX.load(['tOverlap', 'tcc'], [rn, tn])

        coord_transcripts = cgWig.loadSingleWigTranscript(wigDir, chrom, strand, 'transcript')

        for oID in oNX.tOverlap:

                
                tChrom, tStrand, start, end = bioLibCG.tccSplit(oNX.tcc[oID])
                if tStrand == '1':
                        tStrand = '-1'
                else:
                        tStrand = '1'

                if tChrom != chrom or tStrand != strand: continue
                
                oNX.tOverlap[oID] = False
                for i in xrange(start, end + 1):
                        if i in coord_transcripts:
                                oNX.tOverlap[oID] = True
                                break
        

        oNX.save()

def updateContext(oFN, wigDir, chrom, strand, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, cgDegPeak.Peak)
        oNX.load(['tcc', 'context'], [rn, tn])
        
                
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

def updateGSequence(oFN, seqFN, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, cgDegPeak.Peak)
        oNX.load(['gSequence'], [rn, tn])

        f = open(seqFN, 'r')
        i = 0
        for line in f:
                ls = line.strip().split('\t')
                seq = ls[0]
                oNX.gSequence[i] = seq
                i += 1

        oNX.save()

def updateGScore(oFN, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, cgDegPeak.Peak)
        oNX.load(['gScore', 'gSequence'], [rn, tn])

        for oID in oNX.gScore:

                oNX.gScore[oID] = gZipEntropy.gZipEntropy(oNX.gSequence[oID])


        oNX.save()

def updateRepeatStatus(oFN, wigDir, chrom, strand):

        #load oRNAs
        oNX = cgNexusFlat.Nexus(oFN, cgDegPeak.Peak)
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

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
