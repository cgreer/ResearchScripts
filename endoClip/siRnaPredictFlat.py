import bioLibCG as cg
import cgConfig as c
import cgPeaks
import cgGenes
import cgOriginRNAFlat
import cgAlignmentFlat
import compareData as compare
import mapScan
import math
import cgNexusFlat



def updateSequence(oFN, seqFN, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['sequence'], [rn, tn])

        f = open(seqFN, 'r')
        i = 0
        for line in f:
                ls = line.strip().split('\t')
                seq = ls[0]
                oNX.sequence[i] = seq
                i += 1

        oNX.save()

def updateTcc(oFN, tccFN, rn = None, tn = None):

        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['tcc'], [rn, tn])

        f = open(tccFN, 'r')
        i = 0
        for line in f:
                ls = line.strip().split('\t')
                tcc = ls[0]
                oNX.tcc[i] = tcc 
                i += 1

        oNX.save()

def truncateAlignments(alignFN, minNumMismatches, maxNumMismatches, outFN):
        maxNumMismatches = int(maxNumMismatches)
        minNumMismatches = int(minNumMismatches)


        fOut = open(outFN, 'w')
        f = open(alignFN, 'r')
        for line in f:

                #check for mismatches
                numMismatches = int(line.split(' ')[8])
                if numMismatches < minNumMismatches:
                        continue
                
                if numMismatches > maxNumMismatches:
                        continue

                #check that the query is inside the target
                if not (int(line.split(' ')[3]) == int(line.split(' ')[6]) - 1):
                        continue
                if not int(line.split(' ')[2]) == 0:
                        continue

                fOut.write(line)

def updateTargetIDs(oFN, aFN, rn = None, tn = None):
       
        #load the data 
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['targets'], [rn, tn])
        
        #get ids of alignments I need (set)
        oIDs = set()
        for oID in oNX.targets: oIDs.add(oID)

        #load only alignments I need
        c = {'sID' : lambda x: x in oIDs}
        aNX = cgNexusFlat.Nexus(aFN, cgAlignmentFlat.cgAlignment)
        aNX.load(['sID'], conditions = c )

        #clear targets that are there.
        for oID in oNX.targets:
                oNX.targets[oID] = []

        #update the targets for oRNAs
        for aID in aNX.sID:
                oID = aNX.sID[aID]
                oNX.targets[oID].append(aID)
        
        #save
        oNX.save()

def updateTargetIDsFiltered(oFN, aFN, rn = None, tn = None):
        '''CAUTION: NO SELECTION BEING MADE!!!'''

        #load the data 
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['filteredTargets'], [rn, tn])
        
        '''
        #get ids of alignments I need (set)
        oIDs = set()
        for oID in oNX.filteredTargets: oIDs.add(oID)
        '''

        #load only alignments I need
        '''c = {'sID' : lambda x: x in oIDs}'''
        aNX = cgNexusFlat.Nexus(aFN, cgAlignmentFlat.cgAlignment)
        aNX.load(['sID'])

        #clear targets that are there.
        for oID in oNX.filteredTargets:
                oNX.filteredTargets[oID] = []

        #update the targets for oRNAs
        for aID in aNX.sID:
                oID = aNX.sID[aID]
                try:
                        oNX.filteredTargets[oID].append(aID)
                except KeyError: #another process is taking care of this one
                        pass
        
        #save
        oNX.save()

def getEntropy(sequence):

        l = len(sequence)
        #get frequency of every nucleotide
        freqDict = {}
        for i in ['A', 'C', 'T', 'G']:
                freqDict[i] = sequence.count(i)/float(l)

        entropy = 0.0
        for i in ['A', 'C', 'T', 'G']:
                f = freqDict[i]
                if f == 0.0: continue
                entropy += (-1)*f*(math.log(f, 2))

        return entropy

def updateSignificant(resultsFN, simulationAverageFN, outFN):
        
        id_avgNum = {}
        f = open(simulationAverageFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id_avgNum[int(ls[0])] = float(ls[1])
        
        f = open(resultsFN, 'r')
        newLines = []
        for line in f:
                ls = line.strip().split('\t')
                id = int(ls[0])
                numTargets = float(len(ls[4].split(',')))
                try:
                        numExpected = id_avgNum[id]
                except KeyError:
                        numExpected = 0

                sigFlag = 'SIG'
                if numTargets < numExpected:
                        sigFlag = 'NON'

                #Do Calculations here


                updateVal = sigFlag

                #update newLines
                newLines.append(cg.appendToLine(line, updateVal, int(8)))
        f.close()
        
        
        #update file
        f = open(outFN, 'w')
        f.writelines(newLines)
        f.close()


def updateEntropy(oFN, rn = None, tn = None):
	
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['entropy', 'sequence'], [rn, tn])
        
        for oID in oNX.entropy:

                oNX.entropy[oID] = getEntropy(oNX.sequence[oID]) 

        oNX.save()

def updateSmallExpression(oFN, cName, rn = None, tn = None):
	
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['eLevel', 'tcc'], [rn, tn])

        for oID in oNX.eLevel:

	        stretch = cgPeaks.stretch(oNX.tcc[oID], cName) #this stretch contains values for small library...
	        highValue = stretch.getHighestLevel()
	        oNX.eLevel[oID] = highValue

        oNX.save()        

def updateMicroRNAOverlap(aDir, microFN):
	
        oRNA_DC = cgNexusFlat.dataController(aDir, cgOriginRNA.OriginRNA)
	id_oRNA = oRNA_DC.load()
        
        #Put micro and small coords into lists
        microCoords = []
        smallCoords = []
        f = open(microFN, 'r')
        microCoords = [x.strip() for x in f]
        f.close()
        smallCoords = [x.tcc for x in id_oRNA.values()]

        #overlap them
        smallOverlaps = compare.compareTwoTcc(microCoords, smallCoords, 2)


        #For each sRNA, save overlap value.
        for oRNA in id_oRNA.values():
                oRNA.microOverlap = oRNA.tcc in smallOverlaps
	        	
	
        oRNA_DC.commit(id_oRNA)


def updateTargetsExpression(resultsFN, targetsFN, inputPosition, updatePosition, outFN):
	
        #load target expression dict
        f = open(targetsFN, 'r')
        targetsDict = {} # tID: eLevel
        for line in f:
                targetsDict[int(line.strip().split('\t')[0])] = int(line.strip().split('\t')[2])
        f.close()


        #For each sRNA, get target Expression.
	f = open(resultsFN, 'r')
	newLines = []
	for line in f:
		targets = line.strip().split('\t')[int(inputPosition)]
		targets = targets.strip().split(',')
                
                maxExpressionLevel = 0
                totalExpressionLevel = 0
                for tID in targets:
                        tID = int(tID)
                        tExpressionLevel = targetsDict[tID]

                        totalExpressionLevel += targetsDict[tID]
                        if tExpressionLevel > maxExpressionLevel:
                                maxExpressionLevel = tExpressionLevel

	        	
		#update newLines
                newLine = cg.appendToLine(line, maxExpressionLevel, int(updatePosition))
		newLines.append(cg.appendToLine(newLine, totalExpressionLevel, int(updatePosition) + 1))
                
	f.close()
	
	
	#update file
	f = open(outFN, 'w')
	f.writelines(newLines)
	f.close()

def transcriptSetOverlapTargets(aDir):

	geneSetFN = '/home/chrisgre/dataSources/known/Human/geneSets/ensemblAllTranscripts.tsv'
	allExons = cgGenes.createGeneSetFromFile(geneSetFN)

	#get degradome TCCS
	#note that you need to test the AS peaks, this is the location of the targetted transcript
        
        aDC = cgNexusFlat.dataController(aDir, cgAlignment.cgAlignment)
        id_alignment = aDC.load()
        
        #create list of unique tccs.
        uniqTccs = []
        for alignment in id_alignment.values():
                chrom, strand, start, end = cg.tccSplit(alignment.tTcc)
                offset = alignment.tStart
                sLen = alignment.sLength
                if strand == '1':
                        start = start - 19 + offset
                        end = start + sLen
                else:
                        end = end + 19 - offset
                        start = end - sLen

                tcc = cg.makeTcc(chrom, strand, start, end)
                if tcc not in uniqTccs: uniqTccs.append(tcc)

        degTccs = [cg.convertToAS(x) for x in uniqTccs]

	#find all overlapping exons/transcripts, then all results sequences that overlap exons
	overlappingExons = allExons.transcriptOverlaps(degTccs)
        overlappingExonTccs = [x.tcc for x in overlappingExons]
	overlappingDegTccs = compare.compareTwoTcc(degTccs, overlappingExonTccs, 1)

        #update
        for obj in id_alignment.values():         
                chrom, strand, start, end = cg.tccSplit(alignment.tTcc)
                offset = alignment.tStart
                sLen = alignment.sLength

                if strand == '1':
                        start = start - 19 + offset
                        end = start + sLen
                else:
                        end = end + 19 - offset
                        start = end - sLen

                tcc = cg.makeTcc(chrom, strand, start, end)
                degTcc = cg.convertToAS(tcc)

                if degTcc in overlappingDegTccs:
                        obj.transcriptOverlap = True
	        else:
                        obj.transcriptOverlap = False 

        aDC.commit(id_alignment)

def transcriptSetOverlap(aDir, AS):
        AS = bool(AS)

	geneSetFN = '/home/chrisgre/dataSources/known/Human/geneSets/ensemblAllTranscripts.tsv'
	allExons = cgGenes.createGeneSetFromFile(geneSetFN)

	#get degradome TCCS
	#note that you need to test the AS peaks, this is the location of the targetted transcript
        oRNA_DC = cgNexusFlat.dataController(aDir, cgOriginRNA.OriginRNA)
	id_oRNA = oRNA_DC.load()
        if AS == True:
                degTccs = [cg.convertToAS(x.tcc) for x in id_oRNA.values()]
        else:
                degTccs = [x.tcc for x in id_oRNA.values()]

	#find all overlapping exons/transcripts, then all results sequences that overlap exons
	overlappingExons = allExons.transcriptOverlaps(degTccs)
	#print len(overlappingExons), "num of overlapping exons"
        overlappingExonTccs = [x.tcc for x in overlappingExons]
	overlappingDegTccs = compare.compareTwoTcc(degTccs, overlappingExonTccs, 1)


	#write new file
        for obj in id_oRNA.values():         
                if AS:
                        degTcc = cg.convertToAS(obj.tcc)
                else:
                        degTcc = obj.tcc

                if degTcc in overlappingDegTccs:
                        obj.transcriptOverlap = True
	        else:
                        obj.transcriptOverlap = False 

        oRNA_DC.commit(id_oRNA)	
	
def transcriptSetOverlapDegFile(degFile):

	geneSetFN = '/home/chrisgre/dataSources/known/Human/geneSets/ensemblAllTranscripts.tsv'
	allExons = cgGenes.createGeneSetFromFile(geneSetFN)

	#get degradome TCCS
	#note that you need to test the AS peaks, this is the location of the targetted transcript
        
        degTccs = []
        f = open(degFile, 'r')
        for line in f:
                ls = line.strip().split('\t')
                degTccs.append(ls[1])
        f.close()
                        

        degTccs = [cg.convertToAS(x) for x in degTccs]

	#find all overlapping exons/transcripts, then all results sequences that overlap exons
	overlappingExons = allExons.transcriptOverlaps(degTccs)
	#print len(overlappingExons), "num of overlapping exons"
        overlappingExonTccs = [x.tcc for x in overlappingExons]
	overlappingDegTccs = compare.compareTwoTcc(degTccs, overlappingExonTccs, 1)

        
        f = open(degFile, 'r')
	newLines = []
	for line in f:
	        
                degTcc = cg.convertToAS(ls[1])
               
                inTran = '0'
                if degTcc in overlappingDegTccs:
                        inTran = '1'

		#update newLines
                newLine = cg.appendToLine(line, inTran, 3)
                
	f.close()
        
def transcriptSetOverlapDegFileHitmap(degFile, runningChrom, runningStrand):

	geneSetFN = '/home/chrisgre/dataSources/known/Human/geneSets/ensemblAllTranscripts.tsv'
	allExons = cgGenes.createGeneSetFromFile(geneSetFN)
        transcriptTccs = []
        for gene in allExons.set.values():
                for transcript in gene.transcripts:
                        transcriptTccs.append(transcript.tcc)

        #create hitmap
        coordSet = set()
        for tcc in transcriptTccs:
                chrom, strand, start, end = cg.tccSplit(tcc)
                
                if chrom != runningChrom:
                        continue

                if strand != runningStrand:
                        continue

                for i in range(start, end + 1):
                        coordSet.add(i)

        #find overlapping degTccs
        print 'done creating hitmap'
        

        f = open(degFile, 'r')
	newLines = []
	for line in f:
	        ls = line.strip().split('\t') 
                degTcc = cg.convertToAS(ls[1])
                chrom, strand, start, end = cg.tccSplit(degTcc)
                if chrom != runningChrom:
                        continue

                if strand != runningStrand:
                        continue

                inTran = '0'
                for i in xrange(start, end + 1):
                        if i in coordSet:
                                inTran = '1'
                                break

		#update newLines
                newLine = cg.appendToLine(line, inTran, 3)
                newLines.append(newLine)         
	f.close()

        f = open(degFile + '.%s.%s' % (runningChrom, runningStrand), 'w')
        f.writelines(newLines)
        f.close()


if __name__ == "__main__":
	import sys
        cg.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
