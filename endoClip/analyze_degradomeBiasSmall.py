import bioLibCG
import cgNexusFlat
import cgOriginRNAFlat
import cgAlignmentFlat
import cgDegPeak
import cgIContext
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats
import math

def intronBias(degFN, iContextDir, insFN, maxX = 9000):
        maxX = int(maxX)

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
        
        histVals = []
        histVals2 = []
        histVals3 = []
        fOut = open(insFN, 'w')
        for dID in dNX.tcc:

                chrom, strand, start, end = bioLibCG.tccSplit(dNX.tcc[dID])
                dLevel = dNX.eLevel[dID]

                if strand == '1':
                        strand = '-1'
                elif strand == '-1':
                        strand = '1'
                else:
                        print 'WHAT?!'

                for iConID in dNX.iContexts[dID]:

                        if iConID == -1: continue
                        conType = chrom_strand_NX[chrom][strand].type[iConID]
                        if not '3UTR' in conType: continue
                        conTcc = chrom_strand_NX[chrom][strand].tcc[iConID]

                        cChrom, cStrand, cStart, cEnd = bioLibCG.tccSplit(conTcc)
                        iLength = cEnd - cStart

                        if cStrand == '1':
                                dist5 = start - cStart
                                dist3 = cEnd - end
                                pString = [iConID, dID, dNX.tcc[dID], conTcc, dist5, dist3, iLength]
                                fOut.write('\t'.join([str(x) for x in pString]) + '\n')
                        elif cStrand == '-1':
                                dist5 = cEnd - end
                                dist3 = start - cStart
                                pString = [iConID, dID, dNX.tcc[dID], conTcc, dist5, dist3, iLength]
                                fOut.write('\t'.join([str(x) for x in pString]) + '\n')
                        else:
                                print 'WHAT?!'
                        dist5 = [dist5] * dLevel
                        dist3 = [dist3] * dLevel
                        #if dist5[0] < maxX: histVals.extend(dist5)
                        if dist3[0] < maxX: histVals2.extend(dist3)
                        if iLength < maxX: histVals3.append(iLength)
                        
                        #pick only one
                        #break

        #histVals = [math.log(x, 10) if x > 0 else 0 for x in histVals]
        #histVals2 = [math.log(x, 10) if x > 0 else 0 for x in histVals2]
        #histVals3 = [math.log(x, 10) if x > 0 else 0 for x in histVals3]
        plt.xlim(0, maxX)
        plt.hist(histVals, 500, color = 'b', histtype = 'step') 
        plt.hist(histVals2, 500, color = 'g', histtype = 'step') 
        plt.hist(histVals3, 500, color = 'k', histtype = 'step') 
        plt.show()
                        

def intronBiasNormalized(degFN, iContextDir, insFN, minX = 0, maxX = 9000, inType = 'INTRON'):
        maxX = int(maxX)
        minX = int(minX)

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
        
        chrom_strand_ids = {}
        pos_dCount = {}


        conIDs = []
        fOut = open(insFN, 'w')
        countPassed = 0
        for dID in dNX.tcc:

                chrom, strand, start, end = bioLibCG.tccSplit(dNX.tcc[dID])
                dLevel = dNX.eLevel[dID]

                if strand == '1':
                        strand = '-1'
                elif strand == '-1':
                        strand = '1'
                else:
                        print 'WHAT?!'

                for iConID in dNX.iContexts[dID]:

                        if iConID == -1: continue
                        conType = chrom_strand_NX[chrom][strand].type[iConID]
                        if not inType in conType: continue
                        conTcc = chrom_strand_NX[chrom][strand].tcc[iConID]
                        countPassed += 1

                        cChrom, cStrand, cStart, cEnd = bioLibCG.tccSplit(conTcc)
                        iLength = cEnd - cStart
                        if iLength < 110: continue
                        if cStrand == '1':
                                dist5 = start - cStart
                                dist3 = cEnd - end
                                pString = [iConID, dID, dNX.tcc[dID], conTcc, dist5, dist3, iLength]
                                fOut.write('\t'.join([str(x) for x in pString]) + '\n')
                        elif cStrand == '-1':
                                dist5 = cEnd - end
                                dist3 = start - cStart
                                pString = [iConID, dID, dNX.tcc[dID], conTcc, dist5, dist3, iLength]
                                fOut.write('\t'.join([str(x) for x in pString]) + '\n')
                        else:
                                print 'WHAT?!'
                        #if dist5 < maxX: pos_dCount[dist5] = pos_dCount.get(dist5, 0) + 1
                        if minX <= dist3 <= maxX: 
                                pos_dCount[dist3] = pos_dCount.get(dist3, 0) + 1 
                                chrom_strand_ids.setdefault(cChrom, {}).setdefault(cStrand, []).append(iConID)
                                conIDs.append(iConID)              
                                break 
     
        print '# CONIDS in RANGE', len(conIDs)
        for i in range(minX, maxX + 1):
                print i, pos_dCount[i]

        print '---'

        print countPassed
        newCount = 0
        pos_numIntrons = {}
        #make intron length
        for chrom in chrom_strand_ids:
                for strand, idList in chrom_strand_ids[chrom].items():
                        idList = set(idList)
                        newCount += len(idList)
                        for id in idList:
                                introns = chrom_strand_NX[chrom][strand]
                                c, s, st, en = bioLibCG.tccSplit(introns.tcc[id])
                                iLen = en - st
                                for i in xrange(iLen):
                                        if i > maxX: break
                                        pos_numIntrons[i] = pos_numIntrons.get(i, 0) + 1
        print newCount                                        

        #now calculate normalized score for each nt.
        scores = []
        for pos in sorted(pos_dCount.keys()):
                if pos < 0: continue
                dCount = pos_dCount[pos]
                numIntrons = pos_numIntrons[pos]
                score = float(dCount)/float(numIntrons)
                score = dCount
                scores.append(score)
        
        #get num overlapping nts for each CONTEXT
        numPerPos = []
        for pos in sorted(pos_numIntrons.keys()):
                numPerPos.append(pos_numIntrons[pos])

        plt.title(inType)
        plt.subplot(211)
        plt.plot(scores)
        plt.subplot(212)
        plt.plot(numPerPos)

        plt.show()

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
