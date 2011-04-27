import bioLibCG
import compareData
import cgDB
import cgOriginRNA
import cgGenes3
import dumpObj


#Get in terms of tccs
def oneToOne(objectList, objProperty):

        property_object = {}
        for obj in objectList:
                property_object[obj.__dict__[objProperty]] = obj

        return property_object

def oneToMany(objectList, objProperty):

        property_objects = {}
        for obj in objectList:
                property = obj.__dict__[objProperty]
                if property in property_objects:
                        property_objects[property].append(obj)
                else:
                        property_objects[property] = [obj]
        return property_objects

def updateContext(oDir, geneSetFN):
       
        print 'loading oRNA'
        oDC = cgDB.dataController(oDir, cgOriginRNA.OriginRNA)
        id_oRNA = oDC.load()
        print 'loading gene set'
        geneSet = cgGenes3.createGeneSetEditing(geneSetFN)

        
        #Get in terms of tccs
        print 'Joining'
        oTcc_oRNA = oneToOne(id_oRNA.values(), 'tcc')
        tTcc_transcripts = oneToMany(geneSet.transcripts, 'tcc')

        #Overlap tccs
        print 'overlapping'
        oTcc_tTccs = compareData.getIndividualOverlaps(oTcc_oRNA.keys(), tTcc_transcripts.keys(), 1)


        #create final dictionary containing {oRNA : [transcript, ..]}
        oRNA_transcripts = {}
        for oTcc in oTcc_tTccs:
                oRNA = oTcc_oRNA[oTcc]
                oRNA_transcripts[oRNA] = []
                for tTcc in oTcc_tTccs[oTcc]:
                        oRNA_transcripts[oRNA].extend(tTcc_transcripts[tTcc])

        print 'get context info'
        #Go through each site and find out what it overlaps, and if it is in a coding region...
        
        ds = bioLibCG.dominantSpotter(['EXON_INTRON', '3UTR', '5UTR', 'EXON', 'INTRON'])
        for oRNA in oRNA_transcripts:

                oRNA.transcriptIDs = []
                oRNA.transcriptContexts = []
                oRNA.transcriptTypes = []
                oRNA.transcriptCodingTypes = []
                
                if len(oRNA_transcripts[oRNA]) == 0:
                        continue

                for transcript in oRNA_transcripts[oRNA]:
                       
                        codingTranscript =  '_coding' in transcript.tType
                        tType = None
                        codingFlag = None

                        tTypes = [x[1] for x in transcript.getOverlappingElements(oRNA.tcc)]
                        
                        #categorize border types
                       
                        tType = ds.spotItem(tTypes)
                          
                        if tType == 'EXON' or 'EXON_INTRON':
                                if codingTranscript:
                                        codingFlag = 'C'
                                else:
                                        codingFlag = 'NC'
                        else:
                                codingFlag = 'NC'

                        
                        oRNA.transcriptIDs.append(transcript.id)
                        oRNA.transcriptContexts.append(tType)
                        oRNA.transcriptTypes.append(transcript.tType)
                        oRNA.transcriptCodingTypes.append(codingFlag)

        oDC.commit(id_oRNA)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(updateContext, sys.argv)

