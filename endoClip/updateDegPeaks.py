import bioLibCG
import cgDegPeak
import cgAlignmentFlat
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

def updateTOverlapOneRun(oFN, rn = None, tn = None):
        '''Dont need to do it by chromosome because it is small enough'''
        '''Also dont need to flip the strand because the wig is opposite as well'''

        oNX = cgNexusFlat.Nexus(oFN, cgDegPeak.Peak)
        oNX.load(['context', 'tOverlap'], [rn, tn])
        
        for oID in oNX.context:
                
                oNX.tOverlap[oID] = ('INTER' not in oNX.context[oID]) 
       
        oNX.save()

def updateTranscriptOverlap(oFN, wigDir, chrom, strand, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, cgDegPeak.Peak)
        oNX.load(['tOverlap', 'tcc'], [rn, tn])

        #load the AS wig file for this degradome strand
        if strand == '1':
                strand = '-1'
        else:
                strand = '1'
        
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
        
         
        if strand == '1':
                strand = '-1'
        else:
                strand = '1'
        
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

def updateContextAlignment(oFN, wigDir, chrom, strand, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, cgAlignmentFlat.cgAlignment)
        oNX.load(['tTcc', 'context'], [rn, tn])
        
         
        if strand == '1':
                strand = '-1'
        else:
                strand = '1'
        
        print 'loading wig'
        coord_contexts = cgWig.loadSingleWigContext(wigDir, chrom, strand, 'context') 
        print 'done loading'

        ds = bioLibCG.dominantSpotter(['C_EXON', 'C_3UTR', 'C_5UTR', 'NC_EXON', 'NC_3UTR', 'NC_5UTR', 'C_INTRON', 'NC_INTRON', 'INTER']) 


        for oID in oNX.tTcc:

                oChrom, oStrand, start, end = bioLibCG.tccSplit(oNX.tTcc[oID])
                
                #deg wigs is AS to actual clipping site
                if oStrand == '1':
                        oStrand = '-1'
                else:
                        oStrand = '1'
        
                if oChrom == chrom and oStrand == strand:

                        contexts = coord_contexts.get(start, 'INTER').split(',')
                        oNX.context[oID] = ds.spotItem(contexts)

        
        oNX.save()

def updateTypeAlignment(oFN, wigDir, chrom, strand, rn = None, tn = None):
        '''This is for ALIGNMENTS...NOT DEG PEAKS!'''

        oNX = cgNexusFlat.Nexus(oFN, cgAlignmentFlat.cgAlignment)
        oNX.load(['tTcc', 'type'], [rn, tn])
        
                
        if strand == '1':
                strand = '-1'
        else:
                strand = '1'
        
        print 'loading wig'
        coord_types = cgWig.loadSingleWigContext(wigDir, chrom, strand, 'tType') 
        print 'done loading'
        
        domOrder = ['microRNA_noncoding',
        'lincRNA_noncoding',
        'longNC_noncoding',
        'miRNA_pseudogene_noncoding',
        'Mt_rRNA_noncoding',
        'Mt_tRNA_noncoding',
        'Mt_tRNA_pseudogene_noncoding',
        'rRNA_noncoding',
        'rRNA_pseudogene_noncoding',
        'scRNA_pseudogene_noncoding',
        'snoRNA_noncoding',
        'snoRNA_pseudogene_noncoding',
        'snRNA_noncoding',
        'snRNA_pseudogene_noncoding',
        'tRNA_pseudogene_noncoding',
        'pseudogene_noncoding',
        'protein_coding',
        'None']
        
        ds = bioLibCG.dominantSpotter(domOrder)


        for oID in oNX.tTcc:

                oChrom, oStrand, start, end = bioLibCG.tccSplit(oNX.tTcc[oID])
                
                #deg wigs is AS to actual clipping site
                if oStrand == '1':
                        oStrand = '-1'
                else:
                        oStrand = '1'
                
                if oChrom == chrom and oStrand == strand:

                        tranTypes = coord_types.get(start, 'None').split(',')
                        types = [x.split(':')[1] if x != 'None' else 'None' for x in tranTypes]
                        types = list(set(types))
                        oNX.type[oID] = ds.spotItem(types)

        
        oNX.save()

def updateType(oFN, wigDir, chrom, strand, rn = None, tn = None):

        oNX = cgNexusFlat.Nexus(oFN, cgDegPeak.Peak)
        oNX.load(['tcc', 'type'], [rn, tn])
        
                
        if strand == '1':
                strand = '-1'
        else:
                strand = '1'
        
        print 'loading wig'
        coord_types = cgWig.loadSingleWigContext(wigDir, chrom, strand, 'tType') 
        print 'done loading'
        
        domOrder = ['protein_coding', 'miRNA', 'rRNA', 'snoRNA', 'snRNA', 'lincRNA', 'longNC', 'misc_RNA',
                    'Mt_rRNA', 'Mt_tRNA', 'misc_RNA_pseudogene', 'Mt_tRNA_pseudogene', 'polymorphic_pseudogene',
                    'processed_transcript', 'miRNA_pseudogene', 'rRNA_pseudogene', 'scRNA_pseudogene',
                    'snoRNA_pseudogene', 'snRNA_pseudogene', 'tRNA_pseudogene', 'pseudogene', 'TR_C_gene',
                    'TR_J_gene', 'TR_V_gene', 'TR_V_pseudogene', 'IG_C_gene', 'IG_C_pseudogene', 'IG_D_gene',
                    'IG_J_gene', 'IG_J_pseudogene', 'IG_V_gene', 'IG_V_pseudogene', 'unknown', 'None']
        
        ds = bioLibCG.dominantSpotter(domOrder)


        for oID in oNX.tcc:

                oChrom, oStrand, start, end = bioLibCG.tccSplit(oNX.tcc[oID])
                
                #deg wigs is AS to actual clipping site
                if oStrand == '1':
                        oStrand = '-1'
                else:
                        oStrand = '1'
                
                if oChrom == chrom and oStrand == strand:

                        tranTypes = coord_types.get(start, 'None').split(',')
                        print tranTypes
                        types = [x.split(':')[1] if x != 'None' else ['None'] for x in tranTypes]

                        oNX.type = oNX.transcriptType = ds.spotItem(types)

        
        oNX.save()

def updateGSequence(oFN, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, cgDegPeak.Peak)
        oNX.load(['gSequence', 'sequence'], [rn, tn])

        for dID in oNX.sequence:
                oNX.gSequence[dID] = oNX.sequence[dID][16:36]

        oNX.save()

def updateGScore(oFN, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, cgDegPeak.Peak)
        oNX.load(['gScore', 'gSequence'], [rn, tn])

        for oID in oNX.gScore:

                oNX.gScore[oID] = gZipEntropy.gZipEntropy(oNX.gSequence[oID])


        oNX.save()

def updateRepeatCount(oFN, rCountFN, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, cgDegPeak.Peak)
        oNX.load(['repeatCount'], [rn, tn])

        id_count = {}
        f = open(rCountFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id_count[int(ls[0])] = int(ls[1])


        for oID in oNX.repeatCount:
                oNX.repeatCount[oID] = id_count.get(oID, 0)



        oNX.save()

def updateTotalContig(oFN, rn = None, tn = None):

        oNX = cgNexusFlat.Nexus(oFN, cgDegPeak.Peak)
        oNX.load(['sequence', 'totalContig'], [rn, tn])

        for oID in oNX.sequence:
                seq = oNX.sequence[oID] 

                highestLength = 1
                cLength = 1
                letters = list(seq)
                for i,letter in enumerate(letters):
                        if i == 0: continue

                        if letters[i] == letters[i-1]:
                                cLength += 1
                                if cLength > highestLength:
                                        highestLength = cLength
                        else:
                                cLength = 1

                oNX.totalContig[oID] = highestLength

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

def updateIContext(oFN, wigDir, chrom, strand, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, cgDegPeak.Peak)
        oNX.load(['tcc', 'iContexts'], [rn, tn])
        
         
        if strand == '1':
                strand = '-1'
        else:
                strand = '1'
        
        print 'loading wig'
        coord_contexts = cgWig.loadSingleWigContext(wigDir, chrom, strand, 'iContext') 
        print 'done loading'

        for oID in oNX.tcc:

                oChrom, oStrand, start, end = bioLibCG.tccSplit(oNX.tcc[oID])
                
                #deg wigs is AS to actual clipping site
                if oStrand == '1':
                        oStrand = '-1'
                else:
                        oStrand = '1'
        
                if oChrom == chrom and oStrand == strand:

                        contexts = coord_contexts.get(start, '-1').split(',')
                        contexts = [int(x) for x in contexts]

                        oNX.iContexts[oID] = contexts

        
        oNX.save()

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
