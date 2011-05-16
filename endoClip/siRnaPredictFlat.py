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
                numMismatches = int(line.split('\t')[8])
                if numMismatches < minNumMismatches:
                        continue
                
                if numMismatches > maxNumMismatches:
                        continue

                #check that the query is inside the target
                if not (int(line.split('\t')[3]) == int(line.split('\t')[6]) - 1):
                        continue
                if not int(line.split('\t')[2]) == 0:
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

def updateELevel(oFN, wigDir, rn = None, tn = None):
        '''Dont need to do it by chromosome because it is small enough'''
        '''Also dont need to flip the strand because the wig is opposite as well'''

        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['tcc', 'eLevel'], [rn, tn])
        
        wigDict = cgWig.loadWigDict(wigDir)
        
        for oID in oNX.eLevel:
                
                coord_value = cgWig.getExpressionProfile(oNX.tcc[oID], wigDict)
                oNX.eLevel[oID] = max(coord_value.values())

        oNX.save()

def updateContext(oFN, wigDir, chrom, strand, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, cgDegPeak.Peak)
        oNX.load(['tcc', 'context'], [rn, tn])
        
                
        print 'loading wig'
        coord_contexts = cgWig.loadSingleWigContext(wigDir, chrom, strand, 'context') 
        print 'done loading'

        ds = bioLibCG.dominantSpotter(['C_EXON', 'C_3UTR', 'C_5UTR', 'NC_EXON', 'NC_3UTR', 'NC_5UTR', 'C_INTRON', 'NC_INTRON', 'INTER']) 


        for oID in oNX.tcc:

                oChrom, oStrand, start, end = bioLibCG.tccSplit(oNX.tcc[oID])
                if oChrom == chrom and oStrand == strand:

                        contexts = coord_contexts.get(start, 'INTER').split(',')
                        oNX.context[oID] = ds.spotItem(contexts)

        
        oNX.save()


if __name__ == "__main__":
	import sys
        if sys.argv[1] == 'help':
                cg.gd(sys.argv[0])
        else:
                cg.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
