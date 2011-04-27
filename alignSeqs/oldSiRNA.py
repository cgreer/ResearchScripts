

def predictSiRNA(peakFN, cName):
	'''Take peaks and overlap with degradome values.
	degradome values will be in 3' only wig files
	Apparently the degradome data is mapped on wrong strand'''
	
	name = cg.getBaseFileName(peakFN, naked = True)
        pExtend = 6
        cExtend = 3
	f = open(peakFN, 'r')
	smallPeaks = [x.strip() for x in f]
	f.close()
	
	oF = open('%s.degradome' % name, 'w')
	for peakCoord in smallPeaks:
		
                tccChrom, tccStrand, tccStart, tccEnd = cg.tccSplit(peakCoord)

		#get highest degradome values
		stretch = cgPeaks.stretch(peakCoord, cName) #this stretch contains values for degradome...
		highValue = stretch.getHighestLevel()
		highCoord = stretch.getHighestLevel(coord = True)
		
             
                #check if expression in cleaveRange range, get total in that range
                positions = stretch.profile.keys()
                positions.sort()
                cleaveTotal = 0
                

                if tccStrand == '1':
                        for i in range(positions[0] + pExtend + 10 - cExtend, positions[0] + pExtend + 10 + cExtend):# 11 from 5' end 
                                cleaveTotal += stretch.profile[i]
                                #print i, stretch.profile[i]
                else:
                        for i in range(positions[-1] - pExtend - 10 - cExtend, positions[-1] - pExtend - 10 + cExtend):
                                cleaveTotal += stretch.profile[i]
                                #print i, stretch.profile[i]
		
                #output
		oF.write('%s\t.\t%s\t%s\t%s\n' % (peakCoord, highValue, '.', cleaveTotal))

class readInfo:
        def __init__(self, readLine, readType = 'Bowtie'):
               
               if readType == 'Bowtie':
                       self.initBowtie(readLine)

        def initBowtie(self, readLine):
                ls = readLine.strip().split('\t')     
                self.name = ls[0]
                self.strand = cg.switchStrandFormat(ls[1])
                self.chrom = ls[2]
                self.start = int(ls[3])
                self.end = self.start + len(ls[4]) #+1 ?
                self.scores = ls[5]
                if len(ls) > 7:
                        self.mismatches = ls[7]
                else:
                        self.mismatches = None

        def returnMismatchPositions(self):
                
                #parse them first
                mm = []
                if self.mismatches is None:
                        return None
                else:
                        mSplit = self.mismatches.split(',')
                        for mis in mSplit:
                                mm.append(int(mis.split(':')[0]))
                
                pos = []
                #Now calculate positions (absolute)
                if self.strand == '1':
                        for off in mm:
                                pos.append(self.start + off)
                else:
                        for off in mm:
                                pos.append(self.end - off)

                return pos


def markMismatches(inFile, smallDirectory):
        #peak in file
        #look up reads that belong to peak
        #How many have mismatches at position --> put in two columns --> num reads with mis  num reads
        
        f = open(inFile, 'r')
        peakCoords = [line.strip().split('\t')[0] for line in f]
        f.close()
        
        matchDict = {} # peakcoord: [numReads, numMis]
        for peakCoord in peakCoords:
                print peakCoord
                chrom, strand, start, end = cg.tccSplit(peakCoord)
                if chrom == 'chrM': continue
                matchDict[peakCoord] = []
                reads = mapScan.scanCoord(peakCoord, smallDirectory)
                #what if there aren't any reads?
                numReads = len(reads)
                numMisReads = 0
                for read in reads:
                        rInfo = readInfo(read)
                        if rInfo.strand == '1':
                                lowRange = start + 6 + 11 - 3
                                highRange = start + 6 + 11 + 3
                                mPositions = rInfo.returnMismatchPositions()
                                if not mPositions: continue
                                for mm in mPositions:
                                        if lowRange < mm < highRange:
                                                numMisReads += 1
                        else:
                                lowRange = end - 6 - 11 - 3
                                highRange = end - 6 - 11 + 3
                                mPositions = rInfo.returnMismatchPositions()
                                if not mPositions: continue
                                for mm in mPositions:
                                        if lowRange < mm < highRange:
                                                numMisReads += 1
                matchDict[peakCoord].extend([numReads, numMisReads])
        
        #write to new file
        f = open(inFile, 'r')
        newLines = []
        for line in f:
                peakCoord = line.strip().split('\t')[0]
                if peakCoord in matchDict:
                        newLines.append(cg.appendToLine(line, matchDict[peakCoord][0], 6))
                        newLines[-1] = cg.appendToLine(newLines[-1], matchDict[peakCoord][1], 7)
        f.close()

        f = open(inFile + '.mismatch', 'w')
        f.writelines(newLines)
        f.close()


def plotOffset(fN):
	
	f = open(fN, 'r')
	values = [int(x.strip().split('\t')[3]) for x in f if x.strip().split('\t')[3] != 'None']
	#values = [1,2,3,4,5]
	cgPlot.plotHistogram(values) 
