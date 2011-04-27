import bioLibCG
import cgGenes3
import cgEdit
import dumpObj
import cgSeqMod

def updateSynonomous(eFN, gFN, resultsFN, outFN):
                
        #Load Transcripts and Editing Sites
        print 'Loading editing sites'
        eSites = cgEdit.loadEditingSites(eFN)
        print 'Loading gene set'
        geneSet = cgGenes3.createGeneSetEditing(gFN)

        codingTID_eID = {}
        f = open(resultsFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                if ls[4] == 'C':
                        codingTID_eID[ls[2]] = int(ls[0])

        #Get coding Transcripts
        codingTranscripts = {} #tID : eID ! many:one always!
        f = open(resultsFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                if ls[4] == 'C':
                        codingTranscripts[ls[2]] = int(ls[0])


        eID_eSite = {}
        for eSite in eSites:
                eID_eSite[eSite.ID] = eSite

        tID_transcript = {}
        for transcript in geneSet.transcripts:
                tID[transcript.id] = transcript

        codingT_eSite
        for tID in codingTID_eID:
                eID = codingTID_eID[tID]
                t = tID_transcript[tID]
                e = eID_eSite[eID]
                


        print 'Creating scroll dict'
        scrollDict = {} # transcript: eSite
        for tID in codingTranscripts:
                e = eJoint[codingTranscripts[tID]]
                try:
                        t = tJoint[tID]
                        scrollDict[t] = e
                except KeyError:
                        pass

        print 'Deducing synonomous'
        map = cgSeqMod.loadCodonMap('hg19')
        finalDict = {} # tID: [SYN, AAA, AAB, G, A]
        #Figure out if they are synonomous
        for t in scrollDict:

                
                eSite = scrollDict[t]
                #dumpObj.dumpObj(t)
                #dumpObj.dumpObj(eSite)
                
                ePositionInMRNA = t.getRelativePositionMRNA(eSite.coordinate -1)
              
                if ePositionInMRNA == -1:
                        print t.id, 'should not be designated coding...'
                        continue
                
                #grab mRNA and emRNA
                mRNA = t.getMRNA(coding = True)
                emRNA = t.getMRNA(coding = True)

                if mRNA[ePositionInMRNA] != 'A':
                        print 'wrong position', t.id, '%s:%s' % (eSite.chromosome, eSite.coordinate), eSite.strand, mRNA[ePositionInMRNA - 5: ePositionInMRNA -1], mRNA[ePositionInMRNA], mRNA[ePositionInMRNA + 1: ePositionInMRNA + 5]

                #edit the site
                emRNA = list(emRNA)
                emRNA[ePositionInMRNA] = 'G'
                emRNA = ''.join(emRNA)

                #Test the protein sequences
                pRNA = cgSeqMod.translateRNA(mRNA, map)
                epRNA = cgSeqMod.translateRNA(emRNA, map)
                
                #print t.parent, t.id
                newString = ['%s  ' % x for x in list(pRNA)]
                newString = ''.join(newString)

                if pRNA[0] != 'M':
                        print 'Non-canonical Start AA:', pRNA[0:5], mRNA[:10]
                if pRNA[-1] != '*':
                        print 'Non-canonical End AA:', pRNA[-5:], mRNA[-10:]

                

                #compare the codons.

                mCodonList = cgSeqMod.getCodonListFromRNA(mRNA)
                emCodonList = cgSeqMod.getCodonListFromRNA(emRNA)
                compareList = zip(mCodonList, emCodonList)
                synFlag = 'SYN'
                
                codonNumber = ePositionInMRNA // 3
                
                codonPair = compareList[codonNumber]
                print t.id
                print eSite.ID
                print mCodonList[:codonNumber]
                print mRNA[:ePositionInMRNA]
                bCodon = codonPair[0]
                aCodon = codonPair[1]
                
                baa = cgSeqMod.translateRNA(bCodon, map)
                aaa = cgSeqMod.translateRNA(aCodon, map)
                if baa != aaa:
                        synFlag = 'NON'
                        bCodonList = list(bCodon)
                        aCodonList = list(aCodon)
                        matchedLetters = zip(bCodonList, aCodonList)
                        for pair in matchedLetters:
                                if pair[0] != 'A':
                                        if pair[1] == 'G' and pair[0] != 'G':
                                                print 'messed up codon switch', bCodonList, aCodonList
                                                print t.parent, '%s:%s' % (eSite.chromosome, eSite.coordinate), eSite.strand, bCodon, aCodon, baa, aaa
                else:
                        synFlag = 'SYN'

                                

                finalDict[t.id] = [synFlag, bCodon, aCodon, baa, aaa]



        print 'writing to file'
        #update line by line
        newLines = []
        f = open(resultsFN, 'r')
	for line in f:
                newLine = line.strip()
                tID = line.strip().split('\t')[2]
                
                if tID in finalDict:
                        newLine = newLine + '\t%s\t%s\t%s\t%s\t%s\n' % (finalDict[tID][0], finalDict[tID][1],finalDict[tID][2],finalDict[tID][3],finalDict[tID][4])
                else:
                        newLine = newLine + '\tNA\tNA\tNA\tNA\tNA\n'

                newLines.append(newLine)
	f.close()
	
	
	#update file
	f = open(outFN, 'w')
	f.writelines(newLines)
	f.close()

def betterSynonymous(eFN, gFN, contextFN, outFN, refBase = 'A', eBase = 'G'):

        print 'loading e sites'
        eSites = cgEdit.loadEditingSites(eFN)
        print 'loading geneSet'
        geneSet = cgGenes3.createGeneSetEditing(gFN)

        contextInfo = {} # eID: tID : [UTR, C]
        f = open(contextFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                eID = int(ls[0])
                tID = ls[2]
                cInfo = [ls[3], ls[4]]
                if eID not in contextInfo:
                        contextInfo[eID] = {}
                        contextInfo[eID][tID] = cInfo
                else:
                        contextInfo[eID][tID] = cInfo

        eID_tIDs = {}
        f = open(contextFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                eID = int(ls[0])
                tID = ls[2]
                if tID not in eID_tIDs.setdefault(eID, []): eID_tIDs[eID].append(tID)

        eID_eSite = {}
        for eSite in eSites:
                eID_eSite[eSite.ID] = eSite

        tID_transcript = {}
        for transcript in geneSet.transcripts:
                tID_transcript[transcript.id] = transcript

        eSite_transcripts = {}
        for eID in eID_tIDs:
                eSite = eID_eSite[eID]
                tList = []
                for tID in eID_tIDs[eID]:
                        if tID == 'NONE': continue
                        if tID_transcript.get(tID, None) == None: continue
                        tList.append(tID_transcript[tID])
                eSite_transcripts[eSite] = tList

        outF = open(outFN, 'w')
        map = cgSeqMod.loadCodonMap('hg19')
        for eSite in eSite_transcripts:
                
                for transcript in eSite_transcripts[eSite]:
                       
                        siteType, codingType = contextInfo[eSite.ID][transcript.id]
                        if '_noncoding' in transcript.tType:
                                continue
                        if codingType != 'C':
                                continue

                        
                        ePositionInMRNA = transcript.getRelativePositionMRNA(eSite.coordinate -1)
                        mRNA = transcript.getMRNA(coding = True)
                        emRNA = transcript.getMRNA(coding = True)

                        if mRNA[ePositionInMRNA] != refBase:
                                print 'Editing site was not an A...'

                        #edit the site
                        emRNA = list(emRNA)
                        emRNA[ePositionInMRNA] = eBase
                        emRNA = ''.join(emRNA)

                        #Test the protein sequences
                        pRNA = cgSeqMod.translateRNA(mRNA, map)
                        epRNA = cgSeqMod.translateRNA(emRNA, map)
                        if pRNA[0] != 'M':
                                print 'Non-canonical Start AA:', pRNA[0:5], mRNA[:10]
                        if pRNA[-1] != '*':
                                print 'Non-canonical End AA:', pRNA[-5:], mRNA[-10:]

                        #compare the codons.
                        mCodonList = cgSeqMod.getCodonListFromRNA(mRNA)
                        emCodonList = cgSeqMod.getCodonListFromRNA(emRNA)
                        compareList = zip(mCodonList, emCodonList)
                        codonNumber = ePositionInMRNA // 3
                        codonPair = compareList[codonNumber]
                        bCodon = codonPair[0]
                        aCodon = codonPair[1]
                        baa = cgSeqMod.translateRNA(bCodon, map)
                        aaa = cgSeqMod.translateRNA(aCodon, map)
                        synFlag = 'SYN'
                        if baa != aaa:
                                synFlag = 'NON'
                                bCodonList = list(bCodon)
                                aCodonList = list(aCodon)
                                matchedLetters = zip(bCodonList, aCodonList)
                                for pair in matchedLetters:
                                        if pair[0] != 'A':
                                                if pair[1] == 'G' and pair[0] != 'G':
                                                        print 'messed up codon switch', bCodonList, aCodonList
                                                        print t.parent, '%s:%s' % (eSite.chromosome, eSite.coordinate), eSite.strand, bCodon, aCodon, baa, aaa
                        
                                                                
                        outF.write('\t'.join([str(eSite.ID), transcript.parent, transcript.id, synFlag, bCodon, aCodon, baa, aaa]) + '\n')


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(betterSynonymous, sys.argv)
        
