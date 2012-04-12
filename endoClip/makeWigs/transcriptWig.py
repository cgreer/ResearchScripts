import cgWig
import bioLibCG
import compareData

def writeWigDictToWig(wigDict, chrom, strand, assembly, name, outDir, blankValue = 0):
        '''hopefully the coords are in zero based.  And this script will convert it to 0,1'''

        #init
        coords = sorted(wigDict.keys())
        lDict = bioLibCG.returnChromLengthDict(assembly)
        chromEnd = lDict[chrom] 

        outFN = outDir + '/%s.%s.%s.wig' % (name, chrom, strand)
        f = open(outFN, 'w')

        #write first blank line
        f.write('%s\t%s\t%s\t%s\n' % (chrom, 0, coords[0], blankValue))
      
        #p.tell(' 'beginning block', coords[0] 
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
                                #p.tell(' 'writing last equal block', blockStart, prevCoord, prevValue
                                #finish last block, with NO blank block
                                f.write('%s\t%s\t%s\t%s\n' % (chrom, blockStart, prevCoord + 1, prevValue))
                               
                                #init next block
                                prevCoord = coord
                                prevValue = currValue
                                blockStart = coord


                else: #finish last block, write zero block, start another block
                       
                       #last
                       #p.tell(' 'finishing last block', blockStart, prevCoord
                       f.write('%s\t%s\t%s\t%s\n' % (chrom, blockStart, prevCoord + 1, prevValue))
                       
                       #zero
                       #p.tell(' 'zero After:', prevCoord, coord
                       f.write('%s\t%s\t%s\t%s\n' % (chrom, prevCoord + 1, coord, blankValue))

                       #init next block
                       prevCoord = coord
                       blockStart = coord
                       prevValue = currValue
                       
                
        #write last block and last blank block line
        f.write('%s\t%s\t%s\t%s\n' % (chrom, blockStart, coord + 1, prevValue))
        f.write('%s\t%s\t%s\t%s\n' % (chrom, coord + 1, chromEnd, blankValue))
        f.close()

def makeTranscriptWig(tranFN, wigDir, chrom, strand, species = 'hg19'):

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
        writeWigDictToWig(coord_id, chrom, strand, species, 'transcript', wigDir, 'None')       
        
def makeGeneWig(tranFN, wigDir, chrom, strand):

        coord_id = {}
        f = open(tranFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                
                tChrom, tStrand = ls[1], bioLibCG.switchStrandFormat(ls[2])
                if tChrom != chrom or tStrand != strand:
                        continue
                gID = ls[10]
                #gID = gID.replace(" ", "_")
                tStart, tEnd = int(ls[3]), int(ls[4]) - 1

                for i in xrange(tStart, tEnd + 1):
                        coord_id[i] = coord_id.get(i, '')  + '%s$' % gID #$ is used because of spaces

        #unique, string
        for i, ids in coord_id.iteritems():
                coord_id[i] = ','.join([x for x in set(ids.strip().split('$')) if x != ''])

        #write wig to file                
        writeWigDictToWig(coord_id, chrom, strand, 'hg19', 'ALL', wigDir, 'None')       

def makeContextWig(tranFN, wigDir, chrom, strand, species = 'hg19'):

        p = bioLibCG.cgPrint()                               
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
                codingStatus = '_coding' in ls[13]
                tID = ls[0]

                #debug
                p.show = False

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

                #p.tell(tStart, tEnd, cStart, cEnd, exonPairs, intronPairs) 
                
                #take care of messy UTRs and assign utr ranges
                #5UTR
                if strand == '1':
                        if cStart == tStart or cStart == tEnd + 1:
                                p.tell('5 is none')
                                range5 = ()
                        else:
                                range5 = (tStart, cStart - 1)
                else:
                        if cEnd + 1 == tStart or cEnd + 1 == tEnd + 1:
                                p.tell('5 is none')
                                range5 = ()
                        else:
                                range5 = (cEnd + 1, tEnd)

                
                #3UTR
                if strand == '1':
                        if cEnd + 1 == tStart or cEnd + 1 == tEnd + 1:
                                p.tell('3 is none')
                                range3 = ()
                        else:
                                range3 = (cEnd + 1, tEnd)
                else:
                        if cStart == tStart or cStart == tEnd + 1:
                                p.tell('3 is none')
                                range3 = ()
                        else:
                                range3 = (tStart, cStart - 1)
                                
                p.tell('ranges', range5, range3)
                p.tell('intronRange', intronPairs)
                utr5 = compareData.subtractTwoRanges([range5], intronPairs)
                utr3 = compareData.subtractTwoRanges([range3], intronPairs)
                
                p.tell('utr', utr5, utr3)

                p.tell('exon before', exonPairs)
                exonPairs = compareData.subtractTwoRanges(exonPairs, [range5])
                exonPairs = compareData.subtractTwoRanges(exonPairs, [range3])
                p.tell('exon after', exonPairs)

                debugSpot = 23631989 

                #5UTR
                for pair in utr5:
                        p.tell('filling utr5', pair[0], pair[1])
                        for i in xrange(pair[0], pair[1] + 1):
                                if i == debugSpot: p.tell('*** 5UTR', codingStatus, tID)
                                if codingStatus:
                                        coord_id[i] = coord_id.get(i, '')  + 'C_5UTR '
                                else:
                                        coord_id[i] = coord_id.get(i, '')  + 'NC_5UTR '
                
                #Exons
                for pair in exonPairs:
                        p.tell('filling exons', pair[0], pair[1])
                        for i in xrange(pair[0], pair[1] + 1):
                                if i == debugSpot: p.tell('*** exon', codingStatus, tID)
                                if codingStatus:
                                        coord_id[i] = coord_id.get(i, '')  + 'C_EXON '
                                else:
                                        coord_id[i] = coord_id.get(i, '')  + 'NC_EXON '

                #Introns
                for pair in intronPairs:
                        p.tell('filling introns', pair[0], pair[1])
                        for i in xrange(pair[0], pair[1] + 1):
                                if i == debugSpot: p.tell('*** INTRON', codingStatus, tID)
                                if codingStatus:
                                        coord_id[i] = coord_id.get(i, '')  + 'C_INTRON '
                                else:
                                        coord_id[i] = coord_id.get(i, '')  + 'NC_INTRON '

                #3UTR
                for pair in utr3:
                        p.tell('filling utr3', pair[0], pair[1])
                        for i in xrange(pair[0], pair[1] + 1):
                                if i == debugSpot: p.tell(' *** 3UTR', codingStatus, tID)
                                if codingStatus:
                                        coord_id[i] = coord_id.get(i, '')  + 'C_3UTR '
                                else:
                                        coord_id[i] = coord_id.get(i, '')  + 'NC_3UTR '

                p.show = False

        #uniqify, stringify
        for i, ids in coord_id.iteritems():
                coord_id[i] = ','.join([x for x in set(ids.strip().split(' '))])
        
        #p.tell('finalInfo', utr5, exonPairs, utr3)
        #write wig to file                
        writeWigDictToWig(coord_id, chrom, strand, species, 'context', wigDir, 'INTER')       

def makeTypeWig(tranFN, wigDir, chrom, strand, species): 
        '''Using 14th column in transcripts for type info...might want to use something different?'''

        coord_id = {}
        f = open(tranFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                
                tChrom, tStrand = ls[1], bioLibCG.switchStrandFormat(ls[2])
                if tChrom != chrom or tStrand != strand:
                        continue
                tID = ls[0]
                tType = ls[13] 
                tStart, tEnd = int(ls[3]), int(ls[4]) - 1 #0BASE CONVERSION !!! it might have to be 0BASE for making wig...?

                tIDType = '%s:%s' % (tID, tType)       
                for i in xrange(tStart, tEnd + 1):
                        coord_id[i] = coord_id.get(i, '')  + '%s ' % tIDType

        #unique, string
        for i, ids in coord_id.iteritems():
                coord_id[i] = ','.join([x for x in set(ids.strip().split(' '))])

        #write wig to file                
        writeWigDictToWig(coord_id, chrom, strand, species, 'tType', wigDir, 'None')       


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
