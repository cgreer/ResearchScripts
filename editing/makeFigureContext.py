import bioLibCG
import matplotlib.pyplot as plt
from pylab import *
import cgGenes3

def editingSitesPerGene(contextFN):


        f = open(contextFN, 'r')
        gName_eIDs = {}
        for line in f:
                ls = line.strip().split('\t')
                eID = int(ls[0])
                gName = ls[1]
                if gName == 'NONE':
                        continue
                if gName in gName_eIDs:
                        if eID not in gName_eIDs[gName]:
                                gName_eIDs[gName].append(eID)

                else:
                        gName_eIDs[gName] = [eID]
        histVals = []
        for gName in gName_eIDs:
                histVals.append(len(gName_eIDs[gName]))
        plt.title('Editing Sites Per Gene')
        plt.xlabel('Number of Editing Sites')
        plt.ylabel('Number of Genes')
        plt.hist(histVals, 50)
        plt.show()



def makeContextPieBetter(contextFN, gFN, eFN, passedInfo, outFN):
      
        if passedInfo[0] == passedInfo[1]:
                print passedInfo, 'no such thing'
                return 0

        print 'loading geneSet'
        geneSet = cgGenes3.createGeneSetEditing(gFN)
        typeCount = {'EXON': 0, '3UTR': 0, '5UTR': 0, 'INTRON': 0, 'NONG': 0, 'NONT': 0}

        #joint
        tID_transcript = {}
        for transcript in geneSet.transcripts:
                tID_transcript[transcript.id] = transcript

        eID_Info = {}
        f = open(eFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                eID = int(ls[13])
                chrom = ls[0]
                coord = ls[1]
                newCoord = '%s:%s' % (chrom, coord)
                nEdited = ls[4]
                nTotal = ls[5]
                eID_Info[eID] = [newCoord, nEdited, nTotal]

        eID_gName = {}
        eID_tTypes = {}
        f = open(contextFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                eID = int(ls[0])
                gName = ls[1]
                eID_gName[eID] = gName
                tType = ls[3]
                tName = ls[2]
                if tType == 'INTER':
                        continue
                transcript = tID_transcript[tName]
                tInfo = [tType, transcript]

                if eID in eID_tTypes:
                        eID_tTypes[eID].append(tInfo)
                else:
                        eID_tTypes[eID] = [tInfo]

        eID_finalType = {}

        for eID in eID_tTypes:
                highestType = None
                for tInfo in eID_tTypes[eID]:
                        tType = tInfo[0]
                        transcript = tInfo[1]
                        
                        tCoding = True
                        if '_coding' not in transcript.tType:
                                tCoding = False
                        
                        gCoding = True
                        if '_coding' not in transcript.gType:
                                gCoding = False

                        #print transcript.id, tType, tCoding, gCoding
                        
                        if tType == 'EXON':
                                if tCoding:
                                        highestType = 'EXON'
                                        break
                                else:
                                        if gCoding:
                                                if highestType not in ['EXON', 'INTRON', '3UTR', '5UTR']:
                                                        highestType = 'NONT'
                                        else:
                                                if highestType not in ['EXON', 'INTRON', '3UTR', '5UTR', 'NONT']:
                                                        highestType = 'NONG'

                        elif tType == '3UTR':
                                if tCoding:
                                        if gCoding:
                                                highestType = '3UTR'
                                        else:
                                                highestType = '3UTR'
                                else:
                                        if gCoding:
                                                if highestType not in ['EXON', 'INTRON', '3UTR', '5UTR']:     
                                                        highestType = 'NONT'
                                        else:
                                                if highestType not in ['EXON', 'INTRON', '3UTR', '5UTR', 'NONT']:
                                                        highestType = 'NONG'
                        elif tType == '5UTR':
                                if tCoding:
                                        if gCoding:
                                                if highestType not in ['3UTR']:
                                                        highestType = '5UTR'
                                        else:
                                                if highestType not in ['3UTR']:
                                                        highestType = '5UTR'
                                else:
                                        if gCoding:
                                                if highestType not in ['EXON', 'INTRON', '3UTR', '5UTR']:        
                                                        highestType = 'NONT'
                                        else:
                                                if highestType not in ['EXON', 'INTRON', '3UTR', '5UTR', 'NONT']:
                                                        highestType = 'NONG'
                        
                        elif tType == 'INTRON':

                                if tCoding:
                                        if gCoding:
                                                if highestType not in ['3UTR', '5UTR']:
                                                        highestType = 'INTRON'
                                        else:
                                                if highestType not in ['3UTR', '5UTR']:
                                                        highestType = 'INTRON'
                                else:
                                        if gCoding:
                                                if highestType not in ['EXON', 'INTRON', '3UTR', '5UTR']:     
                                                        highestType = 'NONT'
                                        else:
                                                if highestType not in ['EXON', 'INTRON', '3UTR', '5UTR', 'NONT']:
                                                        highestType = 'NONG'
                eID_finalType[eID] = highestType
                typeCount[highestType] += 1


        outF = open(outFN, 'w')
        for eID in eID_finalType:
                outF.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (eID_Info[eID][0].split(':')[0], eID_Info[eID][0].split(':')[1], passedInfo, eID_gName[eID], eID_finalType[eID], eID_Info[eID][1], eID_Info[eID][2]))

        return 0                


        outF = open(outFN, 'w')
        for type in typeCount:
                outF.write('%s\t%s\n' % (type, typeCount[type]))
        return 0 

        #get fractions of each type
        types = ['EXON', '3UTR', '5UTR', 'INTRON', 'NONG', 'NONT']
        fracs = [typeCount['EXON'], typeCount['3UTR'], typeCount['5UTR'], typeCount['INTRON'], typeCount['NONG'], typeCount['NONT']]
        #print fracs
        labels = ['Exons (%s)' % fracs[0], '3\'UTR (%s)' % fracs[1], '5\'UTR (%s)' % fracs[2], 'Introns (%s)' % fracs[3], 'Noncoding Gene (%s)' % fracs[4], 'Noncoding Transcript (%s)' % fracs[5]] 
        theSum = fracs[0] + fracs[1] + fracs[2] + fracs[3] + fracs[4] + fracs[5]
        fracs = [float(x)/theSum for x in fracs]

        #print fracs


        explode=(0.1, 0.1, 0.1, 0.1, .1, .1)
        pie(fracs, explode=explode, labels=labels, autopct='%1.1f%%', shadow=True)
        title('Editing Site Genomic Location', bbox={'facecolor':'1.0', 'pad':10})

        show()

def makeContextPie(contextFN):
        
        #parse values from context Graph
        f = open(contextFN, 'r')
        geneDict = {} # gID: Type
        for line in f:
                ls = line.strip().split('\t')
                gID = ls[0]
                type = ls[2]
                if gID in geneDict:
                        geneDict[gID].append(type)
                else:
                        geneDict[gID] = [type]

        #get fractions of each type
        types = ['INTRON', 'INTER', 'EXON', '3UTR', '5UTR']
        fracs = []
        theSum = 0.0
        for t in types:
                totalNum = 0.0
                for gene in geneDict:
                        lg = len(geneDict[gene])
                        for gType in geneDict[gene]:   
                                if gType == t:
                                        #totalNum += 1.0/lg
                                        #theSum += 1.0/lg
                                        totalNum += 1.0
                                        theSum += 1
                fracs.append(totalNum)


        
        
        labels = ['Introns (%s)' % fracs[0], 'Intergenic(%s)' % fracs[1], 'Exons (%s)' % fracs[2], '3\'UTR (%s)' % fracs[3], '5\'UTR (%s)' % fracs[4]]
        fracs = [float(x)/theSum for x in fracs]

        explode=(0.1, 0.1, 0.2, 0.1, 0.1)
        pie(fracs, explode=explode, labels=labels, autopct='%1.1f%%', shadow=True)
        title('Editing Site Genomic Location', bbox={'facecolor':'1.0', 'pad':10})

        show()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(makeContextPieBetter, sys.argv)

