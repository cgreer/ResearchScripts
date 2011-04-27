import bioLibCG
import cgGenes3
import cgEdit
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

def newExpression(eFN, gFN, contextFN, targetFN, rpkmFN):
        """docstring for newExpression"""

        #print 'loading editing/genes'        
        eSites = cgEdit.loadEditingSites(eFN)
        geneSet = cgGenes3.createGeneSetEditing(gFN)


        #print 'making joints'
        #joints
        eID_eSite = {}
        for eSite in eSites:
                eID_eSite[eSite.ID] = eSite



        gID_eIDs = {}
        eID_tTypes = {}
        f = open(contextFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                gID = ls[1]
                eID = int(ls[0])
                tType = ls[3]

                if gID in gID_eIDs:
                        if eID not in gID_eIDs[gID]:
                                gID_eIDs[gID].append(eID)
                else:
                        gID_eIDs[gID] = [eID]

                if eID in eID_tTypes:
                        eID_tTypes[eID].append(tType)
                else:
                        eID_tTypes[eID] = [tType]

        #print 'updating target sites'
        #update targetting for eSites:
        f = open(targetFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                eID = int(ls[0])
                before = ls[3].split(',')
                if before == ['None']: before = []
                after = ls[4].split(',')
                if after == ['None']: after = []
                eID_eSite[eID].before = before
                eID_eSite[eID].after = after

        #scrolldict
        gene_eSites = {}
        for gID in gID_eIDs:
                try:
                        gene = geneSet.set[gID]
                except KeyError:
                        #print gID, 'not in geneSet'
                        pass

                gene_eSites[gene] = []
                #get eSites
                for eID in gID_eIDs[gID]:
                        eSite = eID_eSite[eID]
                        gene_eSites[gene].append(eSite)

        #print 'updating gene target site info'
        #update before/after editing target sites for GENES
        createdGenes = []
        destroyedGenes = []
        histoVals = []

        for gene in gene_eSites:
                gene.before = []
                gene.after = []
                for eSite in gene_eSites[gene]:
                        if '3UTR' in eID_tTypes[eSite.ID]:
                                for micro in eSite.before:
                                        gene.before.append('%s %s:%s' % (micro, eSite.ID, eSite.coordinate))
                                for micro in eSite.after:
                                        gene.after.append('%s %s:%s' % (micro, eSite.ID, eSite.coordinate))

                if len(gene.before) > len(gene.after):
                        destroyedGenes.append(gene)

                if len(gene.after) > len(gene.before):
                        createdGenes.append(gene)

                change = len(gene.after) - len(gene.before)
                if (len(gene.after) != 0) or (len(gene.before) != 0):
                        histoVals.append(change)

        #print the created/destroyed sites
        for gene in createdGenes:
                print gene.before
                print gene.after

        for gene in destroyedGenes:
                print gene.before
                print gene.after

        '''
        plt.title('Target Site Changes Due To Editing')
        plt.xlabel('Change in Number of Target Sites After Editing')
        plt.ylabel('Number of Genes')
        plt.hist(histoVals, 6)
        plt.show()

        return 0
        '''

        #print out microRNA gene List
        uniqueMicros = {}
        for gene in createdGenes:
                for micro in gene.after:
                        uniqueMicros[micro] = 1
        for gene in destroyedGenes:
                for micro in gene.before:
                        uniqueMicros[micro] = 1

        for micro in uniqueMicros:
                #print micro
                pass

        #print 'loading rpkm'
        gName_ratio = getERatioDict(rpkmFN)
        eChanges = []
        
        for gene in createdGenes:
                try:
                        ratio = gName_ratio[gene.id]
                except KeyError:
                        #print gene.id, 'not in RPKM file --> not expressed'
                        pass
                eChange = math.log(ratio, 2)
                eChanges.append(eChange)
        
        eChanges2 = []
        for gene in destroyedGenes:
                try:
                        ratio = gName_ratio[gene.id]
                except KeyError:
                        #print gene.id, 'not in RPKM file --> not expressed'
                        pass

                eChange = math.log(ratio, 2)
                eChanges2.append(eChange)
        

        #Now plot the histogram
        plt.hist(eChanges, 40, cumulative = True, histtype = 'step', normed = True, label = 'Created')
        plt.hist(eChanges2, 40, cumulative = True, histtype = 'step', normed = True, label = 'Destroyed')
        plt.legend()
        plt.title('eCDF of Genes with Target Site Changes in 3UTRs')
        plt.xlabel('log2(RPKM KD/ RPKM CONTROL)')
        plt.ylabel('Fraction of Genes')
        plt.show()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(newExpression, sys.argv)


