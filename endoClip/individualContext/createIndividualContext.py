import bioLibCG
import cgNexusFlat
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

def makeContextDB(tranFN, chrom, strand, outFN):
        '''outputs cID TYPE TCC, used for making wigs'''

        f = open(tranFN, 'r')
        fOut = open(outFN, 'w')
        id = 0
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
                                range5 = ()
                        else:
                                range5 = (tStart, cStart - 1)
                else:
                        if cEnd + 1 == tStart or cEnd + 1 == tEnd + 1:
                                range5 = ()
                        else:
                                range5 = (cEnd + 1, tEnd)

                
                #3UTR
                if strand == '1':
                        if cEnd + 1 == tStart or cEnd + 1 == tEnd + 1:
                                range3 = ()
                        else:
                                range3 = (cEnd + 1, tEnd)
                else:
                        if cStart == tStart or cStart == tEnd + 1:
                                range3 = ()
                        else:
                                range3 = (tStart, cStart - 1)
                                
                utr5 = compareData.subtractTwoRanges([range5], intronPairs)
                utr3 = compareData.subtractTwoRanges([range3], intronPairs)
                

                exonPairs = compareData.subtractTwoRanges(exonPairs, [range5])
                exonPairs = compareData.subtractTwoRanges(exonPairs, [range3])


                #5UTR
                for pair in utr5:
                        myTcc = bioLibCG.makeTcc(chrom, strand, pair[0] + 1, pair[1] + 1)
                        if codingStatus:
                                pString = [str(id), 'C_5UTR', myTcc]
                                pString = '\t'.join(pString) + '\n'
                                fOut.write(pString)
                        else:
                                pString = [str(id), 'NC_5UTR', myTcc]
                                pString = '\t'.join(pString) + '\n'
                                fOut.write(pString)
                        id += 1                                
                
                #Exons
                for pair in exonPairs:
                        myTcc = bioLibCG.makeTcc(chrom, strand, pair[0] + 1, pair[1] + 1)
                        if codingStatus:
                                pString = [str(id), 'C_EXON', myTcc]
                                pString = '\t'.join(pString) + '\n'
                                fOut.write(pString)
                        else:
                                pString = [str(id), 'NC_EXON', myTcc]
                                pString = '\t'.join(pString) + '\n'
                                fOut.write(pString)
                        id += 1                                

                #Introns
                for pair in intronPairs:
                        myTcc = bioLibCG.makeTcc(chrom, strand, pair[0] + 1, pair[1] + 1)
                        if codingStatus:
                                pString = [str(id), 'C_INTRON', myTcc]
                                pString = '\t'.join(pString) + '\n'
                                fOut.write(pString)
                        else:
                                pString = [str(id), 'NC_INTRON', myTcc]
                                pString = '\t'.join(pString) + '\n'
                                fOut.write(pString)
                        id += 1                                

                #3UTR
                for pair in utr3:
                        myTcc = bioLibCG.makeTcc(chrom, strand, pair[0] + 1, pair[1] + 1)
                        if codingStatus:
                                pString = [str(id), 'C_3UTR', myTcc]
                                pString = '\t'.join(pString) + '\n'
                                fOut.write(pString)
                        else:
                                pString = [str(id), 'NC_3UTR', myTcc]
                                pString = '\t'.join(pString) + '\n'
                                fOut.write(pString)
                        id += 1                                


        f.close()
        fOut.close()

def makeIndividualContextWig(iConFN, wigDir, chrom, strand):
        '''Actual wig file with coordFrom coordTo cID'''

        coord_id = {}
        f = open(iConFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = ls[0]
                tChrom, tStrand, start, end = bioLibCG.tccSplit(ls[2])
                if (chrom != tChrom) or (strand != tStrand): continue 

                for i in xrange(start - 1, end):
                                coord_id[i] = coord_id.get(i, '')  + '%s,' % id


        #uniqify, stringify
        for i, ids in coord_id.iteritems():
                coord_id[i] = ids[:-1] #get rid of trailing comma
        
        #write wig to file                
        writeWigDictToWig(coord_id, chrom, strand, 'hg19', 'iContext', wigDir, 'INTER')       

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
