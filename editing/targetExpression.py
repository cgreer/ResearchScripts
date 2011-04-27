import bioLibCG
import cgEdit
import cgGenes3
import matplotlib.pyplot as plt
import math

def getERatioDict(rpkmFN):
        
        gene_ratio = {}
        f = open(rpkmFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                gName = ls[0]
                eKD = float(ls[1])
                eControl = float(ls[2])

                if (eKD < 1.0) or (eControl < 1.0):
                        continue

                ratio = eKD/eControl
                gene_ratio[gName] = ratio
        f.close()

        return gene_ratio

def makeTargetExpressionHistogram(eFN, targetFN, contextFN, geneFN, eChangeFN):
       

        print 'loading expression ratios'
        gName_eChange = getERatioDict(eChangeFN)

        print 'loading eSites and Transcripts'
        eSites = cgEdit.loadEditingSites(eFN)
        geneSet = cgGenes3.createGeneSetEditing(geneFN)
        
        print 'making joint dicts and loading extra data'
        #joint
        eID_eSite = {}
        for eSite in eSites:
                eID_eSite[eSite.ID] = eSite

        
        #joint
        tID_gName = {}
        for transcript in geneSet.transcripts:
                tID_gName[transcript.id] = transcript.parent


        #load context data
        f = open(contextFN, 'r')
        eID_tID = {} # eID: tID
        tID_tType = {}
        for line in f:
                ls = line.strip().split('\t')
                eID = int(ls[0])
                tID = ls[1]
                tType = ls[2]
                tID_tType[tID] = tType
                if eID in eID_tID:
                        eID_tID[eID].append(tID)
                else:
                        eID_tID[eID] = [tID]
        f.close()
        

        print 'analyzing'
        #Get created or destroyed
        f = open(targetFN, 'r')
        altered = []
        for line in f:
                ls = line.strip().split('\t')
                #created/destroyed
                if (ls[3] != 'None') and (ls[4] == 'None'):
                        altered.append(int(ls[0]))
        f.close()
        
                                       
        print 'number of created/destroyed sites:', len(altered)

        alteredSites = []
        for id in altered:
                alteredSites.append(eID_eSite[id])


        eChanges = []
        gDone = []
        #Get gene names for each eSite
        for eSite in alteredSites:
                genes = []
                for tID in eID_tID[eSite.ID]:
                        if tID == 'NONE':
                                continue

                        gName = tID_gName[tID]
                        if tID_tType[tID] != '3UTR':
                                print 'Not 3UTR', tID
                                continue
                        if gName not in genes:
                                genes.append(gName)

                
                if len(genes) > 1:
                        print 'more than one gene for eSite...', genes
                        continue
                
                if gName in gDone:
                        continue
                else:
                        gDone.append(gName)


                #Now add expression to HistoGram List...
                if gName in gName_eChange:
                        eChange = gName_eChange[gName]
                else:
                        print 'gene not in expression list', gName
                        continue
                eChange = gName_eChange[gName]
                eChange = math.log(eChange, 2)
                eChanges.append(eChange)

     
        #Now plot the histogram
        plt.hist(eChanges, 40)
        plt.xlabel('log2(RPKM KD/ RPKM CONTROL)')
        plt.ylabel('# Genes')
        plt.show()


        
        

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(makeTargetExpressionHistogram, sys.argv)






                  




