import bioLibCG
import cgNexusFlat
import cgOriginRNAFlat
import cgAlignmentFlat
import cgDegPeak
import cgIContext
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats
import scipy
import math
from cgAutoCast import autocast
import random
import GenomeFetch as gf

def updateHitMap(chrom_strand_coord, tcc):

    chrom, strand, start, end = bioLibCG.tccSplit(tcc)

    for i in xrange(start, end + 1):
        chrom_strand_coord.setdefault(chrom, {}).setdefault(strand, set()).add(i)

def checkHitOverlap(chrom_strand_coord, tcc):

    chrom, strand, start, end = bioLibCG.tccSplit(tcc)

    if chrom not in chrom_strand_coord:
        return False

    if strand not in chrom_strand_coord[chrom]:
        return False

    for i in xrange(start, end + 1):
        if i in chrom_strand_coord[chrom][strand]:
            return True

    return False        
        

@autocast
def intronBiasNormalized(degFN, occupancyFN, iContextDir, insFN, minX = 0, maxX = 9000, inType = 'INTRON', prime5 = True, switchStrand = True, removeDups = True, removeMulti = False):

        #load degNX
        dNX = cgNexusFlat.Nexus(degFN, cgDegPeak.Peak)
        dNX.load(['tcc', 'iContexts', 'eLevel'])

        #load all context NXs
        print 'loading NXs'
        chrom_strand_NX = {}
        for chrom in bioLibCG.humanChromosomes:
                for strand in ('1', '-1'):
                        
                        iFN = iContextDir + '/%s.%s.ids' % (chrom, strand) 
                        NX = cgNexusFlat.Nexus(iFN, cgIContext.IContext)
                        NX.load(['type', 'tcc'])

                        chrom_strand_NX.setdefault(chrom, {})[strand] = NX
        print 'done...'
        
        chrom_strand_ids = {}
        pos_dCount = {}
        cLengths = []
        cProportions = []
        conIDs = []
        fOut = open(insFN, 'w')
        countPassed = 0
        hitChrom_strand_coord = {}
        overlapCount = 0
        multiCount = 0
        pickCount = 0
        usedConTccs = set()
        for dID in dNX.ids:

                peakChrom, peakStrand, peakStart, peakEnd = bioLibCG.tccSplit(dNX.tcc[dID])
                dLevel = dNX.eLevel[dID]

                if removeDups:
                    #check for overlaping peaks...
                    if checkHitOverlap(hitChrom_strand_coord, dNX.tcc[dID]):
                        overlapCount += 1
                        continue
                    else:
                        updateHitMap(hitChrom_strand_coord, dNX.tcc[dID]) 


                #hela data is backwards
                if switchStrand:
                    if peakStrand == '1':
                            peakStrand = '-1'
                    elif peakStrand == '-1':
                            peakStrand = '1'
                    else:
                            print 'WHAT?!'


                if len(dNX.iContexts[dID]) == 1:
                    if dNX.iContexts[dID][0] == -1:
                        continue
                
                if removeMulti:  
                    dTypes = [chrom_strand_NX[peakChrom][peakStrand].type[x] for x in dNX.iContexts[dID]]
                    multPass = [inType in cType for cType in dTypes]
                    if multPass.count(True) > 1:
                        multiCount += 1
                        continue

                #check all the context ids this peak overlaps (good use of random Nexus)
                overlappingIDs = dNX.iContexts[dID]
                random.shuffle(overlappingIDs)
                for iConID in overlappingIDs:
                        
                        #get info/filter
                        if iConID == -1: continue
                        conType = chrom_strand_NX[peakChrom][peakStrand].type[iConID]
                        if not inType in conType: continue
                        conTcc = chrom_strand_NX[peakChrom][peakStrand].tcc[iConID]
                        countPassed += 1

                        #get 3/5 end distance
                        cChrom, cStrand, cStart, cEnd = bioLibCG.tccSplit(conTcc)
                        iLength = cEnd - cStart
                        if iLength < 50: continue
                        if cStrand == '1':
                                dist5 = peakStart - cStart
                                dist3 = cEnd - peakEnd
                                pString = [iConID, dID, dNX.tcc[dID], conTcc, dist5, dist3, iLength]
                                fOut.write('\t'.join([str(x) for x in pString]) + '\n')
                        elif cStrand == '-1':
                                dist5 = cEnd - peakEnd
                                dist3 = peakStart - cStart
                                pString = [iConID, dID, dNX.tcc[dID], conTcc, dist5, dist3, iLength]
                                fOut.write('\t'.join([str(x) for x in pString]) + '\n')
                        else:
                                print 'WHAT?!'
                        
                        if conTcc not in usedConTccs:
                            propStat = dist3/float(iLength)
                            if 25.8 <= dist3 <= 26.2:
                                #print '%s:%s-%s %s | %s:%s-%s %s | peakID %s conID %s' % (cChrom, cStart, cEnd, cStrand, peakChrom, peakStart, peakEnd, peakStrand, dID, iConID )
                                print bioLibCG.tccToGB(dNX.tcc[dID]), peakStrand, conTcc
                                pickCount += 1
                            cProportions.append(dist3/float(iLength))
                            usedConTccs.add(conTcc)


                        #add to peak count at that position 
                        if prime5:
                            if minX <= dist5 <= maxX: 
                                    pos_dCount[dist5] = pos_dCount.get(dist5, 0) + 1 
                                    #chrom_strand_ids.setdefault(cChrom, {}).setdefault(cStrand, []).append(iConID)
                                    conIDs.append(iConID)              
                                    cLengths.append(iLength)
                                    break 
                        else:
                            if minX <= dist3 <= maxX: 

                                    pos_dCount[dist3] = pos_dCount.get(dist3, 0) + 1 
                                    #chrom_strand_ids.setdefault(cChrom, {}).setdefault(cStrand, []).append(iConID)
                                    conIDs.append(iConID)              
                                    cLengths.append(iLength)
                                    break
                                     
     
        print '# CONIDS in RANGE', len(conIDs)
        print 'overlap count', overlapCount
        print 'multicount', multiCount
        print 'pickcount', pickCount
        print '---'

        if True:        
            #get occupancy for positions
            pos_numIntrons = {}
            f = open(occupancyFN, 'r')
            for line in f:
                ls = line.strip().split('\t')
                pos, numContexts = int(ls[0]), int(ls[1])
                if pos > maxX: break 
                pos_numIntrons[pos] = numContexts
            f.close()

            #now calculate normalized score for each nt.
            scores = []
            for pos in range(maxX):
                    if pos < 0: continue
                    dCount = pos_dCount.get(pos, 0)
                    numIntrons = pos_numIntrons[pos]
                    score = float(dCount)/float(numIntrons)
                    scores.append(score)
            
            #get num overlapping nts for each CONTEXT
            numPerPos = []
            for pos in sorted(pos_numIntrons.keys()):
                    numPerPos.append(pos_numIntrons[pos])

            primeText = ''
            if prime5:
                primeText = '5P'
            else:
                primeText = '3P'

            #dem plots 
            plt.subplot(411)
            plt.title('%s BIAS (%s)\nlength median/std %s/%s' % (primeText, inType, scipy.median(cLengths), scipy.std(cLengths)) )
            plt.xlabel('# peaks/# introns occupancy')
            plt.plot(scores)
            plt.xlim([0, 1000])

            plt.subplot(412)
            plt.xlabel('# context @ position')
            plt.plot(numPerPos)

            plt.subplot(413)
            plt.xlabel('# context with Y length')
            plt.hist(cLengths, bins = 2000)
            plt.xlim([0,20000])

        n_count = {}
        for n in cProportions:
            n_count[n] = n_count.get(n, 0) + 1

        #for n in sorted(n_count.keys()):
            #if -.001 < n < .02:
                #print n, n_count[n]

        plt.subplot(414)
        plt.xlabel('# context with Y percentage from edge')
        plt.xlim([-.1, 1.1])
        plt.hist(cProportions, bins = 4000)

        print scipy.median(cLengths), scipy.std(cLengths)
        plt.show()

@autocast
def getNumberContextPosition(iContextDir, outFN, cType = 'INTRON', prime5 = True, ):

        #load all context NXs
        print 'loading NXs'
        chrom_strand_NX = {}
        for chrom in bioLibCG.humanChromosomes:
                for strand in ('1', '-1'):
                        iFN = iContextDir + '/%s.%s.ids' % (chrom, strand) 
                        NX = cgNexusFlat.Nexus(iFN, cgIContext.IContext)
                        NX.load(['type', 'tcc'])

                        chrom_strand_NX.setdefault(chrom, {})[strand] = NX
        print 'done...'


        print 'getting context occupancy'
        #get the amount of introns/etc that occupy X
        pos_numIntrons = {}
        for chrom in chrom_strand_NX:
            for strand in chrom_strand_NX[chrom]:
                NX = chrom_strand_NX[chrom][strand]
                for id in NX.ids:
                    if cType not in NX.type[id]: continue 
                    c, s, st, en = bioLibCG.tccSplit(NX.tcc[id])
                    iLen = en - st
                    for i in xrange(0, iLen):
                            pos_numIntrons[i] = pos_numIntrons.get(i, 0) + 1

        fOut = open(outFN, 'w')
        for pos, numIntrons in pos_numIntrons.iteritems():
            fOut.write('%s\t%s\n' % (pos, numIntrons) )
        fOut.close()

def getSeqEnrichment(seqs):

    nt_count = {}
    for seq in seqs:
        for nt in ('A', 'T', 'C', 'G'):
            nt_count[nt] = nt_count.get(nt, 0) + seq.count(nt)

    
    return nt_count 

def getBiasedSeqs(fN, assembly, switchStrand = True):

    seqs = []
    f = open(fN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        chrom, strand, start, end = bioLibCG.tccSplit(ls[0])
        if switchStrand:
            strand = bioLibCG.switchStrand(strand)
        start -= 10
        end += 10 
        seqs.append(bioLibCG.makeTcc(chrom,strand,start,end))
    f.close()

    myG = gf.GenomeFetch(assembly)
    sequences = []
    for i, seq in enumerate(seqs):
        sequences.append(myG.getSequence(seq))
        print '>blah_%s' % i 
        print sequences[-1]

    for let, count in getSeqEnrichment(sequences).items():
        print let, count 

@autocast
def collectIntronSeqs(fN, outFN, assembly, amount = 100, prime3 = True):

    myG = gf.GenomeFetch(assembly)
    fOut = open(outFN, 'w')
    f = open(fN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        if not 'C_INTRON' in ls[1]: continue
        chrom, strand, start, end = bioLibCG.tccSplit(ls[2])

        if end - start < amount: continue
        
        if strand == '1':
            if prime3:
                tcc = bioLibCG.makeTcc(chrom, strand, end - amount, end)
            seq = myG.getSequence(tcc)
            fOut.write('%s\n' % seq[::-1])

        else:
            if prime3:
                tcc = bioLibCG.makeTcc(chrom, strand, start, start + amount)
            seq = myG.getSequence(tcc)
            fOut.write('%s\n' % seq[::-1])


    f.close()
    fOut.close()

@autocast
def plotInfoNTFrameEnrichment(dir, frameWidth, outFN):

    frameNum_seqs = {}
    for fChrom in bioLibCG.humanChromosomes:
        for fStrand in ('1', '-1'):
            print fChrom, fStrand
            fN = '%s/%s.%s.primeSeqs' % (dir, fChrom, fStrand)
            f = open(fN, 'r')
            for line in f:
                ls = line.strip().split('\t')
                seq = ls[0]
                for i, frame in enumerate(bioLibCG.returnFrames(seq, frameWidth)):
                    frameNum_seqs.setdefault(i, []).append(frame)
            f.close()


    let_frameCounts = {}
    for fNum in sorted(frameNum_seqs.keys()):
        seqs = frameNum_seqs[fNum]
        let_count = getSeqEnrichment(seqs)
        for let, count in let_count.items():
            let_frameCounts.setdefault(let, []).append(count)


    fOut = open(outFN, 'w')
    for let, fCounts in let_frameCounts.items():
        fCounts = ','.join([str(x) for x in fCounts])
        fOut.write('%s\t%s\n' % (let, fCounts))
    fOut.close()

def plotNTFrameEnrichment(fN, frameWidth):    
        
    let_frameCounts = {}
    f = open(fN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        let = ls[0]
        counts = [int(x) for x in ls[1].split(',')]
        let_frameCounts[let] = counts
    f.close()
        
    for let, fCounts in let_frameCounts.items():
        plt.plot(fCounts, label = let)

    plt.legend()
    plt.xlabel('frameNumber(%s)' % frameWidth)
    plt.ylabel('# of NT')
    plt.title('3P Intron NT Frame Counts')
    plt.show()

        

           
    


if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
