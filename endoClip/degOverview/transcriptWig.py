import cgWig
import bioLibCG
import compareData

def writeWigDictToWig(wigDict, chrom, strand, assembly, name, outDir, blankValue = 0):
        '''hopefully the coords are in zero based.  And this script will convert it to 0,1'''

        #init
        coords = sorted(wigDict.keys())
        lDict = bioLibCG.returnChromLengthDict('hg19')
        chromEnd = lDict[chrom] 

        outFN = outDir + '/%s.%s.%s.wig' % (name, chrom, strand)
        f = open(outFN, 'w')

        #write first blank line
        f.write('%s\t%s\t%s\t%s\n' % (chrom, 0, coords[0], blankValue))
      
        #print 'beginning block', coords[0] 
        prevCoord = coords[0]
        prevValue = wigDict[coords[0]]
        blockStart = prevCoord
        coords = coords[1:]
        for coord in coords:
                currValue = wigDict[coord]
     
                if coord - 1 == prevCoord:
                        #Does the value differ?
                        if currValue == prevValue:
                                #keep extending block
                                prevCoord = coord
                                prevValue = currValue
                        else:
                                #print 'writing last equal block', blockStart, prevCoord, prevValue
                                #finish last block, with NO blank block
                                f.write('%s\t%s\t%s\t%s\n' % (chrom, blockStart, prevCoord + 1, prevValue))
                               
                                #init next block
                                prevCoord = coord
                                prevValue = currValue
                                blockStart = coord


                else: #finish last block, write zero block, start another block
                       
                       #last
                       #print 'finishing last block', blockStart, prevCoord
                       f.write('%s\t%s\t%s\t%s\n' % (chrom, blockStart, prevCoord + 1, prevValue))
                       
                       #zero
                       #print 'zero After:', prevCoord, coord
                       f.write('%s\t%s\t%s\t%s\n' % (chrom, prevCoord + 1, coord, blankValue))

                       #init next block
                       prevCoord = coord
                       blockStart = coord
                       prevValue = currValue
                       
                
        #write last block and last blank block line
        f.write('%s\t%s\t%s\t%s\n' % (chrom, blockStart, coord + 1, prevValue))
        f.write('%s\t%s\t%s\t%s\n' % (chrom, coord + 1, chromEnd, blankValue))
        f.close()

def makeTranscriptWig(tranFN, wigDir, chrom, strand):

        coord_id = {}
        f = open(tranFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                
                tChrom, tStrand = ls[1], bioLibCG.switchStrandFormat(ls[2])
                if tChrom != chrom or tStrand != strand:
                        continue
                tID = ls[0]
                tStart, tEnd = int(ls[3]), int(ls[4]) - 1

                for i in xrange(tStart, tEnd + 1):
                        coord_id[i] = coord_id.get(i, '')  + '%s ' % tID

        #unique, string
        for i, ids in coord_id.iteritems():
                coord_id[i] = ','.join([x for x in set(ids.strip().split(' '))])

        #write wig to file                
        writeWigDictToWig(coord_id, chrom, strand, 'hg19', 'transcript', wigDir, 'None')       
        
def makeContextWig(tranFN, wigDir, chrom, strand):

        coord_id = {}
        f = open(tranFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                
                tChrom, tStrand = ls[1], bioLibCG.switchStrandFormat(ls[2])
                if tChrom != chrom or tStrand != strand:
                        continue
                tStart, tEnd = int(ls[3]), int(ls[4]) - 1
                cStart, cEnd = int(ls[5]), int(ls[6]) - 1
                exonStarts = [int(x) for x in ls[8][:-1].split(',')]
                exonEnds = [int(x) - 1 for x in ls[9][:-1].split(',')]
                exonPairs = zip(exonStarts, exonEnds)
                codingStatus = '_coding' in ls[15]
                tID = ls[0]
                intronPairs = []
                i = 0
                for pair in exonPairs:
                        if i == 0:
                                i += 1
                                continue
                        iStart = exonPairs[i -1][1] + 1
                        iEnd = exonPairs[i][0] - 1
                        intronPairs.append((iStart, iEnd))
                        i += 1

                
                #take care of messy UTRs and assign utr ranges
                if cStart == tStart or cStart == tEnd + 1:
                        #print '5 is none'
                        range5 = ()
                else:
                        if strand == '1':
                                range5 = (tStart, cStart - 1)
                        else:
                                range5 = (cEnd + 1, tEnd)
                
                if cEnd + 1 == tStart or cEnd + 1 == tEnd + 1:
                        #print '3 is none'
                        range3 = ()
                else:
                        if strand == 1:
                                range3 = (cEnd + 1, tEnd)
                        else:
                                range3 = (tStart, cStart - 1)
                
                #print 'ranges', range5, range3
                #print 'intronRange', intronPairs
                utr5 = compareData.subtractTwoRanges([range5], intronPairs)
                utr3 = compareData.subtractTwoRanges([range3], intronPairs)
                
                #print 'utr', utr5, utr3

                #print 'exon before', exonPairs
                exonPairs = compareData.subtractTwoRanges(exonPairs, [range5])
                exonPairs = compareData.subtractTwoRanges(exonPairs, [range3])
                #print 'exon after', exonPairs

         
                #print 'info', utr5, exonPairs, utr3

                #5UTR
                for pair in utr5:
                        #print 'filling utr5', pair[0], pair[1]
                        for i in xrange(pair[0], pair[1] + 1):
                                if codingStatus:
                                        coord_id[i] = coord_id.get(i, '')  + 'C_5UTR '
                                else:
                                        coord_id[i] = coord_id.get(i, '')  + 'NC_5UTR '
                
                #Exons
                for pair in exonPairs:
                        #print 'filling exons', pair[0], pair[1]
                        for i in xrange(pair[0], pair[1] + 1):
                                if codingStatus:
                                        coord_id[i] = coord_id.get(i, '')  + 'C_EXON '
                                else:
                                        coord_id[i] = coord_id.get(i, '')  + 'NC_EXON '

                #Introns
                for pair in intronPairs:
                        #print 'filling introns', pair[0], pair[1]
                        for i in xrange(pair[0], pair[1] + 1):
                                if codingStatus:
                                        coord_id[i] = coord_id.get(i, '')  + 'C_INTRON '
                                else:
                                        coord_id[i] = coord_id.get(i, '')  + 'NC_INTRON '

                #3UTR
                for pair in utr3:
                        #print 'filling utr3', pair[0], pair[1]
                        for i in xrange(pair[0], pair[1] + 1):
                                if codingStatus:
                                        coord_id[i] = coord_id.get(i, '')  + 'C_3UTR '
                                else:
                                        coord_id[i] = coord_id.get(i, '')  + 'NC_3UTR '

        #uniqify, stringify
        for i, ids in coord_id.iteritems():
                coord_id[i] = ','.join([x for x in set(ids.strip().split(' '))])
        
        #print 'finalInfo', utr5, exonPairs, utr3
        #write wig to file                
        writeWigDictToWig(coord_id, chrom, strand, 'hg19', 'context', wigDir, 'INTER')       

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
