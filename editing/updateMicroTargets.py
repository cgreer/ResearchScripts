import bioLibCG
import cgEdit
import cgMicroRNA
import GenomeFetch
import cgSeqMod
import dumpObj
import cgGenes3


def allindices(string, sub, listindex=[], offset=0):
        i = string.find(sub, offset)
        while i >= 0:
                listindex.append(i)
                i = string.find(sub, i + 1)
        return listindex

def updateLocationBasedTargets(editFN, contextFN, miLocationFN, gFN):

        eSites = cgEdit.loadEditingSites(editFN)
        cgEdit.updateContextEditingSites(eSites, contextFN) #puts the UTR, EXON in eSite.context
        
        geneSet = cgGenes3.createGeneSetEditing(gFN)

        tName_t = {}
        for t in geneSet.transcripts:
                tName_t[t.id] = t


        tName_miInfo = {}
        f = open(miLocationFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                tName = ls[0]
                miName = ls[1]
                loc = int(ls[2])
                tName_miInfo.setdefault(tName, []).append([miName, loc])


        for eSite in eSites:
                if '3UTR' not in eSite.context:
                        continue
                for tName in eSite.transcripts:
                        if tName in tName_miInfo:
                                
                                t = tName_t[tName]
                                
                                for info in tName_miInfo[tName]:
                                        
                                        
                                        miName = info[0]
                                        loc = info[1]

                                        #get the position of e site in mrna for this transcript
                                        ePosition = t.getRelativePositionMRNA(eSite.coordinate, coding = False)
                                        
                                        print tName, miName, loc, ePosition
                                        if loc - 22 <= ePosition <= loc:
                                                print tName, miName, '%s:%s' % (eSite.chromosome, eSite.coordinate)
                                                pass

def checkIfSeedPresent(gName_micros, gFN):

        gf = GenomeFetch.GenomeFetch('hg19')

        print 'loading gene set'
        geneSet = cgGenes3.createGeneSetEditing(gFN)
        print '....done loading'

        outF = open('utrSeeds', 'w')
        done = {}
        for gName in gName_micros:
                micros = gName_micros[gName]
                for transcript in geneSet.set[gName].transcripts:
                        
                        utrCoords = []
                        if len(transcript.utr3) == 0:
                                continue
                        for utrPair in transcript.utr3:
                                utrCoords.extend(utrPair)
                        utrCoords.sort()
                        start, end = utrCoords[0], utrCoords[1]
                        chrom = transcript.chromosome
                        strand = transcript.strand

                        checkSeq = gf.get_seq_from_to(chrom, start, end, strand)
                        print ''
                        print transcript.id
                        print checkSeq
                        checkSeq = checkSeq.replace('T', 'U')

                        for micro in micros:
                                
                                findings = checkSeq.find(micro.comSeed)
                                if findings != -1:
                                        outF.write('%s\t%s\t%s\n' % ( transcript.id, micro.name, transcript.parent))
                                        found = findings
                                        
def get3UTRSeq(transcript):
        gf = GenomeFetch.GenomeFetch('hg19')       
        utrCoords = []
        if len(transcript.utr3) == 0:
                print 'NO UTR'
                return 0
        for utrPair in transcript.utr3:
                utrCoords.extend(utrPair)
        utrCoords.sort()


        starts = [(utrCoords[i*2])for i in range(len(utrCoords)//2)] 
        ends = [utrCoords[1 + i*2] for i in range(len(utrCoords)//2)]
        chrom = transcript.chromosome
        strand = transcript.strand
        if strand == '-1':
                utrCoords.reverse()
        
        for pair in zip(starts, ends):
                start, end = pair[0], pair[1]

                checkSeq = checkSeq + gf.get_seq_from_to(chrom, start, end, strand)
        checkSeq = checkSeq.replace('T', 'U')

        return checkSeq

def checkSeeds(editFN, contextFN, miLocationFN, miSequenceFN, gFN):

        eSites = cgEdit.loadEditingSites(editFN)
        cgEdit.updateContextEditingSites(eSites, contextFN) #puts the UTR, EXON in eSite.context
        
        geneSet = cgGenes3.createGeneSetEditing(gFN)

        tName_t = {}
        for t in geneSet.transcripts:
                tName_t[t.id] = t

        
        miName_miSequence = {}
        f = open(miSequenceFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                name = ls[0]
                seq = ls[1]
                name = 'hsa-' + name
                miName_miSequence[name] = seq

        tName_miInfo = {}
        f = open(miLocationFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                tName = ls[0]
                miName = ls[1]
                loc = int(ls[2])
                tName_miInfo.setdefault(tName, []).append([miName, loc])

        foundIt = []
        notFoundIt = []
        for tName in tName_miInfo:
                
                try:
                        t = tName_t[tName]
                except:
                        continue
                checkSeq = get3UTRSeq(t)
                try:
                        mRNA = t.getMRNA()
                except:
                        continue
                for miInfo in tName_miInfo[tName]:
                
                        miName = miInfo[0]
                        loc = miInfo[1]
                        try:
                                miSequence = miName_miSequence[miName]
                                miSeed = miSequence[1:8]
                        except:
                                continue

                        rcMiSeed = cgSeqMod.reverseComplementSequence(miSeed, True)
                        
                        newLoc = loc - (len(mRNA) - len(checkSeq)) 
                        finding = checkSeq.find(rcMiSeed, newLoc - 25)    
                        if finding != -1:
                                if (0 < newLoc - finding < 30):
                                        newResult = '%s\t%s\t%s\t%s\t%s' % (miName, tName, finding, newLoc, loc)
                                        if newResult not in foundIt: foundIt.append(newResult)
                        else:
                                        
                                        if miName == 'hsa-miR-21':
                                                print loc, len(checkSeq), len(mRNA)
                                                print mRNA
                                                print checkSeq
                                               
                                        newResult = '%s\t%s\t%s\t%s\t%s' % (miName, tName, finding, newLoc, loc)
                                        if newResult not in notFoundIt: notFoundIt.append(newResult)

        print len(foundIt)
        print len(notFoundIt)
        print ''
        for i in foundIt:
                print i
        print ''
        for i in notFoundIt:
                print i

def updateValidatedMicroTargets(editFN, microTargetFN, microSequenceFN, outFN, gFN):
        flankAmount = 6
        eSites = cgEdit.loadEditingSites(editFN)
        cgEdit.updateContextEditingSites(eSites)

        miNames_micros = cgMicroRNA.loadMicroRNAFromValidated(microTargetFN, microSequenceFN)
        gf = GenomeFetch.GenomeFetch('hg19')
        
        #update flanking region
        for eSite in eSites:

                chrom = eSite.chromosome
                coord = eSite.coordinate
                strand = eSite.strand

                flankingSeq = gf.get_seq_from_to(chrom, coord - flankAmount, coord + flankAmount, strand)
                eFlankingSeq = flankingSeq[:flankAmount] + 'G' + flankingSeq[flankAmount + 1:]

                eSite.flank = flankingSeq.replace('T', 'U')
                eSite.eFlank = eFlankingSeq.replace('T', 'U')

        #joint
        gName_micros = {}
        for micro in miNames_micros.values():
                for target in micro.targetGenes:
                        if micro not in gName_micros.setdefault(target, []): gName_micros[target].append(micro)

        gene_m = {}       
        for eSite in eSites:

                sharedMicros = gName_micros.get(eSite.gene)
                if sharedMicros is None:
                        continue
                for micro in sharedMicros:
                        print '' 
                        print micro.name, eSite.gene, micro.sequence, micro.seed
                        print micro.comSeed
                        print eSite.flank, eSite.eFlank
                        print '%s:%s' % (eSite.chromosome, eSite.coordinate), eSite.flank, eSite.gene
                        if micro.comSeed == None:
                                #dumpObj.dumpObj(micro)
                                #print 'miR not in sequence file:', micro.name
                                #print micro.targetGenes
                                continue

                        if eSite.gene in gene_m:
                                if micro not in gene_m[eSite.gene]:
                                        gene_m[eSite.gene].append(micro)
                        else:
                                gene_m[eSite.gene] = [micro]

                        #flanking
                        if micro.comSeed in eSite.flank:
                                eSite.before.append(micro.name)

                        if micro.comSeed in eSite.eFlank:
                                eSite.after.append(micro.name)

        print len(gene_m)
        count = 0
        for g in gene_m:
                print g
                for m in gene_m[g]:
                        print '...', m.name
                        count += 1
        print count
        
        #check if these seeds are in the 
        checkIfSeedPresent(gene_m, gFN) 
                
     
        #write contents to file...
        outF = open(outFN, 'w')
        for eSite in eSites:
                if len(eSite.before) == 0:
                        targets = 'None'
                else:
                        targets = ','.join(eSite.before)

                if len(eSite.after) == 0:
                        eTargets = 'None'
                else:
                        eTargets = ','.join(eSite.after)

                outF.write('%s\t%s:%s\t%s\t%s\t%s\n' % (eSite.ID, eSite.chromosome, eSite.coordinate, eSite.strand, targets, eTargets))
        

def countBySeed(editFN, microFN, flankAmount, outFN):
        '''FlankF amount should be +/- 6'''
        flankAmount = int(flankAmount)

        eSites = cgEdit.loadEditingSites(editFN)
        micros = cgMicroRNA.loadMicroRNAFromTargetScan(microFN, 'hsa')
        gf = GenomeFetch.GenomeFetch('hg19')

        for eSite in eSites:

                chrom = eSite.chromosome
                coord = eSite.coordinate
                strand = eSite.strand

                flankingSeq = gf.get_seq_from_to(chrom, coord - flankAmount, coord + flankAmount, strand)
                eFlankingSeq = flankingSeq[:flankAmount] + 'G' + flankingSeq[flankAmount + 1:]
                flankingSeq.replace('T', 'U')
                eFlankingSeq.replace('T','U')

                checkID = 'hsa-miR-330-5p'
                
                for microRNA in micros:
                                
                        comSeed = cgSeqMod.reverseComplementSequence(microRNA.seed, True)

                        if comSeed in flankingSeq:
                                eSite.microTargets.append(microRNA.id)
                                if microRNA.id == checkID:
                                        print '%s:%s' % (eSite.chromosome, eSite.coordinate), eSite.strand, eSite.gene
                                microRNA.numBefore += 1
                                #print '@', 'flank', eSite.ID, microRNA.id, flankingSeq, comSeed                                

                        if comSeed in eFlankingSeq:
                                microRNA.numAfter += 1
                                if microRNA.id == checkID:
                                        print '%s:%s' % (eSite.chromosome, eSite.coordinate), eSite.strand, eSite.gene
                                eSite.eMicroTargets.append(microRNA.id)
                                #print '@', 'eFlank', eSite.ID, microRNA.id, eFlankingSeq, comSeed
                                
        for micro in micros:
                if micro.numBefore > 0 or micro.numAfter > 0:
                        #print micro.id, micro.numBefore, micro.numAfter
                        pass

        #write contents to file...
        outF = open(outFN, 'w')
        for eSite in eSites:
                if len(eSite.microTargets) == 0:
                        targets = 'None'
                else:
                        targets = ','.join(eSite.microTargets)

                if len(eSite.eMicroTargets) == 0:
                        eTargets = 'None'
                else:
                        eTargets = ','.join(eSite.eMicroTargets)

                outF.write('%s\t%s:%s\t%s\t%s\t%s\n' % (eSite.ID, eSite.chromosome, eSite.coordinate, eSite.strand, targets, eTargets))


def updateMicroTargets(editFN, microFN, flankAmount, outFN):
        '''FlankF amount should be +/- 6'''
        flankAmount = int(flankAmount)

        eSites = cgEdit.loadEditingSites(editFN)
        micros = cgMicroRNA.loadMicroRNAFromTargetScan(microFN, 'hsa')
        gf = GenomeFetch.GenomeFetch('hg19')

        for eSite in eSites:

                chrom = eSite.chromosome
                coord = eSite.coordinate
                strand = eSite.strand

                flankingSeq = gf.get_seq_from_to(chrom, coord - flankAmount, coord + flankAmount, strand)
                eFlankingSeq = flankingSeq[:flankAmount] + 'G' + flankingSeq[flankAmount + 1:]
                

                for microRNA in micros:
                        
                        comSeed = cgSeqMod.reverseComplementSequence(microRNA.seed, True)

                        if comSeed in flankingSeq:
                                eSite.microTargets.append(microRNA.id)
                                print '@', 'flank', eSite.ID, microRNA.id, flankingSeq, comSeed                                

                        if comSeed in eFlankingSeq:
                                eSite.eMicroTargets.append(microRNA.id)
                                print '@', 'eFlank', eSite.ID, microRNA.id, eFlankingSeq, comSeed
                                

                
        #write contents to file...
        outF = open(outFN, 'w')
        for eSite in eSites:
                if len(eSite.microTargets) == 0:
                        targets = 'None'
                else:
                        targets = ','.join(eSite.microTargets)

                if len(eSite.eMicroTargets) == 0:
                        eTargets = 'None'
                else:
                        eTargets = ','.join(eSite.eMicroTargets)

                outF.write('%s\t%s:%s\t%s\t%s\t%s\n' % (eSite.ID, eSite.chromosome, eSite.coordinate, eSite.strand, targets, eTargets))


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(countBySeed, sys.argv)



