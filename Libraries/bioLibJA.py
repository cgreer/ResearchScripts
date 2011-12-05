import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast
import random
import time
import os 


'''Prints n random lines from the file.'''
def randLines(fN, numLines):
    numLines = int(numLines)
    f = open(fN, 'r')
    fileSize = os.path.getsize(fN)
    for x in range(numLines):
        r = random.randint(0,fileSize)
        f.seek(r,0)
        character = f.read(1)
        while character != '\n':
            f.seek(-2, 1)
            character = f.read(1)
        print f.readline()
    f.close()

def randLinesToFile(fN, numLines):
    numLines = int(numLines)
    f = open(fN, 'r')
    fileSize = os.path.getsize(fN)
    for x in range(numLines):
        r = random.randint(0,fileSize)
        f.seek(r,0)
        character = f.read(1)
        while character != '\n':
            f.seek(-2, 1)
            character = f.read(1)
        print f.readline(),
    f.close()

'''From blat, takes the coordinate relative to the read and converts it into genomic coordinate, given the alignment start coordinates and alignment length. Returns the genomic coordinate.''' 
def readToGenomicCoord(alignStarts, alignLen, readCoord):
    readCoord = int(readCoord)
    alignStarts = [int(x) for x in alignStarts]

    alignLen = [int(x) for x in alignLen]
    array = zip(alignStarts, alignLen)
    alignEndGen = [sum(x) for x in array]
    #print alignEndGen

    juncStart = [sum(alignLen[:n+1]) for n in range(len(alignLen))]
    #print juncStart
    array2 = zip(alignStarts[1:], alignEndGen[:-1])
    #print 'array2=', array2

    juncLen = [(x-y) for x,y in array2]
    #print juncLen
    
    totalIntrons = 0
    for i, start in enumerate(juncStart):
        #print i, start
        if readCoord > start:
            # print 'True'
            totalIntrons = totalIntrons + juncLen[i]
    genCoord = alignStarts[0] + totalIntrons + readCoord

    return genCoord
    #print genCoord


'''Takes the mismatch string (MD:Z:mm) from SAM format and the read sequence and returns a list of the mismatch coordinate relative to the read, the reference base (Before), and the read base (After).'''
def parseMismatchStr(mmString, readSeq):
    mmString = mmString.split(':')
    mmStr = mmString[2]        

    After = []
    mmReadCoord = []
    Before = []
    numVar = ''

    for character in mmStr: 
        if character.isdigit() == True:
            numVar = numVar + character
        elif character.isalpha() == True:
            After.append(character.upper())
            numVar = int(numVar)
            mmReadCoord.append(numVar)
            numVar = ''
        else:
            raise NameError(("Something wrong in the syntax '%s'") % mismatch)
    
    mmReadCoord = [i + 1 for i in mmReadCoord]
    mmReadCoord = [sum(mmReadCoord[:n + 1]) for n in range(len(mmReadCoord))]

    for number in mmReadCoord:
        Before.append(readSeq[number - 1])    #0-based index		
    #print 'mmreadcoord', mmReadCoord
    return [mmReadCoord, Before, After]


'''Builds a reference dictionary of all the SNPs to which the mismatches cross-reference. chrom:coord:before:after'''
def buildSNPdict(fSNP):
    print 'Building SNP Dictionary'
    f = open(fSNP, 'r')
    
    chrom_coord_before_after = {}

    for line in f:
        ls = line.strip().split(' ')
        #print ls
        chrom = ls[0]
        coord = int(ls[1])
        ss = ls[2].split('>')
        before = ss[0]
        after = ss[1]
        
        #chrom_coord_before_after.setdefault(chrom, {}).setdefault(coord, {}).setdefault(before, []).append(after)
        chrom_coord_before_after.setdefault(chrom, {}).setdefault(coord, tuple([before, after]))

    f.close()
    return chrom_coord_before_after


'''Corrects mismatch calls that are actually SNPs in the blat alignment. Outputs to a new file. fIN is the uncorrected blat parsed file.'''
def SNPcorrectionBlat(fIN, correctOUT, correctedIDsOUT, fSNP):    
    #print 'Starting SNP correction'
    f = open(fIN, 'r')
    q = open(correctOUT, 'w')
    r = open(correctedIDsOUT, 'w')
    newLines = []
    chrom_coord_before_after = buildSNPdict(fSNP)
    print 'Starting SNP correction'
    #stats
    correctionCount = 0
    noCorrection = 0
    corrected = []
    #pdb.set_trace();
    for num, line in enumerate(f):
        ls = line.strip().split('\t')    
        chrom = ls[5]
        startPos = ls[3]
        numMM = int(ls[20])

        if numMM == 0:
            newLines.append('%s\t%s\t%s\n' % ('\t'.join(ls[0:21]), '0', '\t'.join(ls[21:])))
            noCorrection += 1
            continue
        else:
            mmCoord = ls[21]
            mmCoordSplit = [int(x) for x in mmCoord.strip().split(',')]
            mmBase = ls[22]
            mmBaseSplit = mmBase.strip().split(',')
            readBase = ls[23]
            readBaseSplit = readBase.strip().split(',')
            #mmqScore = ls[18]
            #mmqScoreSplit = mmqScore.strip().split(',')

        #check if lengths of mismatch info are equal.
        a = len(mmCoordSplit) == len(mmBaseSplit) == len(readBaseSplit)

        if a == False:
            print line
            print 'Error: Lists do not match in length.'
            break
        else: 
            'Something is wrong with your lists...'
        
        #SNPannotation = []
        for i in range(len(mmCoordSplit)):
            a = mmCoordSplit[i]
            b = mmBaseSplit[i]
            c = readBaseSplit[i]
            try:
                x = (b,c) == chrom_coord_before_after[chrom][a]
                if x:
                    print line
                    corrected.append('%s\n' %(ls[0]))
                    correctionCount += 1
                    #SNPannotation.append('%s\t%s\t%s\t%s\n' %(('\t'.join(ls)), a, b, c))
                    mmCoordSplit[i] = '!'
                    mmBaseSplit[i] = '!'
                    readBaseSplit[i] = '!'
                    #numMM -= 1
                if x == False:
                    noCorrection += 1
            except KeyError: 
                noCorrection += 1
                pass

        mmCoordSplit = [x for x in mmCoordSplit if x != '!']    
        mmBaseSplit = [x for x in mmBaseSplit if x != '!']
        readBaseSplit = [x for x in readBaseSplit if x != '!']

        numMM2 = len(mmCoordSplit)    
        mmCoordSplit = [str(x) for x in mmCoordSplit]
        #SNPannotation = [str(x) for x in SNPannotation]
        mmCoordSplit = ','.join(mmCoordSplit)
        mmBaseSplit = ','.join(mmBaseSplit)
        readBaseSplit = ','.join(readBaseSplit)
        #mmqScoreSplit = ','.join(mmqScoreSplit)

        newLines.append('%s\t%s\t%s\t%s\t%s\n' %(('\t'.join(ls[:-3])), numMM2, mmCoordSplit, mmBaseSplit, readBaseSplit))  
    q.writelines(newLines)
    #r.writelines(corrected)

    print "Number of SNP corrections = ", correctionCount        
    print "Number of alignments with no corrections = ", noCorrection
    print "Total = ", (correctionCount + noCorrection)
    f.close()
    q.close()


'''Corrects mismatch calls that are actually SNPs in the bowtie.tran alignment. Outputs to a new file. fIN is the uncorrected blat parsed file.'''
def SNPcorrectionBowtieTran(fIN, correctOUT, correctedIDsOUT, fSNP):    
    #print 'Starting SNP correction'
    f = open(fIN, 'r')
    q = open(correctOUT, 'w')
    r = open(correctedIDsOUT, 'w')
    newLines = []
    chrom_coord_before_after = buildSNPdict(fSNP)
    print 'Starting SNP correction'

    #stats
    correctionCount = 0
    noCorrection = 0
    corrected = []

    for num, line in enumerate(f):
        ls = line.strip().split('\t')    
        tcc = ls[3]
        tccsp = tcc.strip().split(':')
        chrom = tccsp[2]
        numMM = int(ls[11])

        if numMM == 0:
            newLines.append('%s\t%s\t%s\n' % ('\t'.join(ls[0:12]), '0', '\t'.join(ls[12:])))
            noCorrection += 1
            continue
        else:
            mmCoord = ls[12]
            mmCoordSplit = [int(x) for x in mmCoord.strip().split(',')]
            mmBase = ls[13]
            mmBaseSplit = mmBase.strip().split(',')
            readBase = ls[14]
            readBaseSplit = readBase.strip().split(',')
            #mmqScore = ls[18]
            #mmqScoreSplit = mmqScore.strip().split(',')

        #check if lengths of mismatch info are equal.
        a = len(mmCoordSplit) == len(mmBaseSplit) == len(readBaseSplit)

        if a == False:
            print line
            print 'Error: Lists do not match in length.'
            break
        else: 
            'Something is wrong with your lists...'
        
        #SNPannotation = []
        for i in range(len(mmCoordSplit)):
            a = mmCoordSplit[i]
            b = mmBaseSplit[i]
            c = readBaseSplit[i]
            try:
                x = (b,c) == chrom_coord_before_after[chrom][a]
                if x:
                    print line
                    corrected.append('%s\n' %(ls[0]))
                    correctionCount += 1
                    #SNPannotation.append('%s\t%s\t%s\t%s\n' %(('\t'.join(ls)), a, b, c))
                    mmCoordSplit[i] = '!'
                    mmBaseSplit[i] = '!'
                    readBaseSplit[i] = '!'
                    #numMM -= 1
                if x == False:
                    noCorrection += 1
            except KeyError: 
                noCorrection += 1
                pass

        mmCoordSplit = [x for x in mmCoordSplit if x != '!']    
        mmBaseSplit = [x for x in mmBaseSplit if x != '!']
        readBaseSplit = [x for x in readBaseSplit if x != '!']

        numMM2 = len(mmCoordSplit)    
        mmCoordSplit = [str(x) for x in mmCoordSplit]
        #SNPannotation = [str(x) for x in SNPannotation]
        mmCoordSplit = ','.join(mmCoordSplit)
        mmBaseSplit = ','.join(mmBaseSplit)
        readBaseSplit = ','.join(readBaseSplit)
        #mmqScoreSplit = ','.join(mmqScoreSplit)

        newLines.append('%s\t%s\t%s\t%s\t%s\n' %(('\t'.join(ls[:-3])), numMM2, mmCoordSplit, mmBaseSplit, readBaseSplit))  
    q.writelines(newLines)
    #r.writelines(corrected)

    print "Number of SNP corrections = ", correctionCount        
    print "Number of alignments with no corrections = ", noCorrection
    print "Total = ", (correctionCount + noCorrection)
    f.close()
    q.close()


''' Removes mapped reads with more than 12 mismatches for RNA editing'''
@autocast
def removeMoreThan12mm(fIN, fOUT, numMMindex):
    print 'Removing reads with more than 12 mismatches...'
    f = open(fIN, 'r')
    q = open(fOUT, 'w')
    newLines = []
    removed = 0

    for line in f:
        ls = line.strip().split('\t')
        numMM = int(ls[numMMindex])
        if numMM > 12:
            removed += 1
            continue
        else:
            newLines.append(line)

    q.writelines(newLines)
    f.close()
    q.close()
    print 'Removed %d reads with more than 12 mismatches.' % removed


if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
