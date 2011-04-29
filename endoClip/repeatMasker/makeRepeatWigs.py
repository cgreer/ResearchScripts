import bioLibCG

def writeSetToWig(wigSet, chrom, strand, assembly, name, outDir):


        print 'if TP in set', (208148750 in wigSet)

        #init
        coords = sorted(wigSet)
        lDict = bioLibCG.returnChromLengthDict('hg19')
        chromEnd = lDict[chrom] 

        outFN = outDir + '/%s.%s.%s.wig' % (name, chrom, strand)
        f = open(outFN, 'w')

        #write first 0 line
        f.write('%s\t%s\t%s\t%s\n' % (chrom, 0, coords[0], 0))
      
        
        prevCoord = coords[0]
        blockStart = prevCoord
        coords = coords[1:]
        for coord in coords:
                if coord - 1 == prevCoord:
                        #keep extending block
                        prevCoord = coord
                else: #finish last block, write zero block, start another block
                       
                       #last
                       f.write('%s\t%s\t%s\t%s\n' % (chrom, blockStart, prevCoord + 1, 1))
                       
                       #zero
                       f.write('%s\t%s\t%s\t%s\n' % (chrom, prevCoord + 1, coord, 0))

                       #init next block
                       prevCoord = coord
                       blockStart = coord
                       
                
        #write last block and last 0 block line
        f.write('%s\t%s\t%s\t%s\n' % (chrom, blockStart, coord + 1, 1))
        f.write('%s\t%s\t%s\t%s\n' % (chrom, coord + 1, chromEnd, 0))
        f.close()



def makeRepeatWigs(selectedChrom, selectedStrand, wigFN, outDir):
        '''meant to be ran in parallel'''
        '''The repeats are strandless'''

        #get repeat locations
        rLocations = set()
        f = open(wigFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                chrom, strand = ls[5], ls[9]
                if chrom != selectedChrom:
                        continue
                start, end = int(ls[6]), int(ls[7])
                
                for i in range(start, end + 1):
                        rLocations.add(i)
                

        
        name = 'REPEAT'
        writeSetToWig(rLocations, selectedChrom, selectedStrand, 'hg19', name, outDir)       

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
