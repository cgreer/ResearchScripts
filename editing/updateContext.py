import bioLibCG 
import compareData
import cgGenes3
import cgEdit
import dumpObj


def updateContext(editFN, geneSetFN, outFN, refBase = 'A'):
        print refBase
        #Load Transcripts and Editing Sites
        print 'loading editing sites'
        eSites = cgEdit.loadEditingSites(editFN, refBase)
        print 'loading gene set'
        geneSet = cgGenes3.createGeneSetEditing(geneSetFN)

        
        #make the eSites 0 based
        for eSite in eSites:
                #redo coordinate and tcc
                eSite.coordinate = eSite.coordinate - 1
                eSite.tcc = bioLibCG.makeTcc(eSite.chromosome, eSite.strand, eSite.coordinate, eSite.coordinate)

        #Create Joint dictionaries
        print 'creating joint dictionaries'
        eJoint = {} #tcc : eSite
        for eSite in eSites:
                eJoint[eSite.tcc] = eSite

        tJoint = {} # tcc : [transcript, ...]
        for transcript in geneSet.transcripts:
                if transcript.tcc in tJoint:
                        tJoint[transcript.tcc].append(transcript)
                else:
                        tJoint[transcript.tcc] = [transcript]

        #Overlap tccs
        print 'overlapping joints'
        ##make new 0-based keys
        tccOverlaps = compareData.getIndividualOverlaps(eJoint.keys(), tJoint.keys(), 1)


        print 'creating final dictionary'
        #create final dictionary containing {edit sites : [transcript, ..]}
        eSiteTranscripts = {} # edit site: [transcript, ..]
        for eTcc in tccOverlaps:
                eSite = eJoint[eTcc]
                eSiteTranscripts[eSite] = []
                for tTcc in tccOverlaps[eTcc]:
                        eSiteTranscripts[eSite].extend(tJoint[tTcc])

        print 'get context info'
        #Go through each site and find out what it overlaps, and if it is in a coding region...
        fOut = open(outFN, 'w')
        for eSite in eSiteTranscripts:
                if len(eSiteTranscripts[eSite]) == 0:
                        #label intergenic
                        tType = 'INTER'
                        codingFlag = 'NC'
                        fOut.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (eSite.ID, 'NONE', 'NONE', tType, codingFlag, 'NONE'))
                        continue

                for transcript in eSiteTranscripts[eSite]:
                       
                        codingTranscript =  '_coding' in transcript.tType
                        tType = None
                        codingFlag = None
                        tTypes = [ x[1] for x in transcript.getOverlappingElements(eSite.tcc)]
                        

                        if '3UTR' in tTypes:
                                tType = '3UTR'
                        elif '5UTR' in tTypes:
                                tType = '5UTR'
                        else:
                                tType = tTypes[0] #has to be one thing...exon or intron
                        #This only works because UTR takes precedence over EXON in TYPE.
                        if tType == 'EXON':
                                if codingTranscript:
                                        codingFlag = 'C'
                                else:
                                        codingFlag = 'NC'
                        else:
                                codingFlag = 'NC'


                        fOut.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (eSite.ID, transcript.parent, transcript.id, tType, codingFlag, transcript.tType))
                        #fOut.write('%s:%s:%s\t%s\t%s\t%s\t%s\n' % (eSite.chromosome, eSite.strand, eSite.coordinate, transcript.parent, transcript.id, tType, codingFlag))

        fOut.close()


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(updateContext, sys.argv)

