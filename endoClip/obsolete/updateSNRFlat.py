import bioLibCG
import cgOriginRNAFlat
import os
import cgNexusFlat

def updateAvgNumTargets(oFN):

        bn = os.path.basename(oFN)
        print 'basename', bn

        oID_numTargets = {}

        for i in range(0,10):
                
                #simFN = '/home/chrisgre/scripts/simulations/simsk50FilteredMasked/simulation.%s/%s' % (i, bn)
                simFN = '/home/chrisgre/scripts/simulations/simsk50/simulation.%s/%s' % (i, bn)
                print simFN
                osNX = cgNexusFlat.Nexus(simFN, cgOriginRNAFlat.OriginRNA)
                osNX.load(['filteredTargets'])
                for oID in osNX.filteredTargets:
                        currTargets = oID_numTargets.get(oID, 0)
                        oID_numTargets[oID] = currTargets + len(osNX.filteredTargets[oID])
        
        #now save it
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['avgNumSimulationTargets'])

        for oID in oNX.avgNumSimulationTargets:
                totalNum = oID_numTargets.get(oID, 0)
                avgNum = float(totalNum)/float(10.0)
                oNX.avgNumSimulationTargets[oID] = avgNum

        oNX.save()

def updateAvgNumSS(oFN):

        bn = os.path.basename(oFN)
        print 'basename', bn

        oID_numSS = {}

        numSims = 100

        print 'getting avg for %s simulations' % numSims
        for i in range(0,numSims):
                
                #simFN = '/home/chrisgre/scripts/simulations/simsk50FilteredMasked/simulation.%s/%s' % (i, bn)
                #simFN = '/home/chrisgre/scripts/simulations/simsk50Fix/simulation.%s/%s' % (i, bn)
                #simFN = '/home/chrisgre/scripts/simulations/mm9/simulation.%s/%s' % (i, bn)
                #simFN = '/home/chrisgre/scripts/simulations/hg19.hela/simulation.%s/%s' % (i, bn)
                simFN = '/home/chrisgre/scripts/simulations/hg19.U87/simulation.%s/%s' % (i, bn)
                osNX = cgNexusFlat.Nexus(simFN, cgOriginRNAFlat.OriginRNA)
                osNX.load(['numSignificantSequences'])
                for oID in osNX.numSignificantSequences:
                        oID_numSS[oID] = oID_numSS.get(oID, 0) + osNX.numSignificantSequences[oID]
        
        #now save it
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['avgNumSS'])

        for oID in oNX.avgNumSS:
                totalNum = oID_numSS.get(oID, 0)
                avgNum = float(totalNum)/float(numSims)
                oNX.avgNumSS[oID] = avgNum
        oNX.save()

def countForError(oDir, filteredFile):

        fList = []
        f = open(filteredFile, 'r')
        for line in f:
                ls = line.strip().split('\t')
                fList.append(int(line.strip()))
       
        oID_numTargets = {}

        for i in range(0,10):
                print i
                simDirRNA = '/home/chrisgre/scripts/simulations/simsk50Filtered/simulation.%s/oRNA' % i
                oDC = cgNexusFlat.dataController(simDirRNA, cgOriginRNA.OriginRNA)
                id_sRNA = oDC.load()
                
                groupTotal = 0
                for id, sRNA in id_sRNA.items():
                        if not id in fList:
                                continue
                        groupTotal +=  len(sRNA.filteredTargets)
                print groupTotal                        
        
def updateSNR(oFN):
        
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['avgNumSimulationTargets', 'filteredTargets', 'snr'])

        for oID in oNX.snr:
                actualNum = len(oNX.filteredTargets[oID])
                avgNum = oNX.avgNumSimulationTargets[oID]
                
                if avgNum == 0: avgNum = 1
                
                SNR = float(actualNum)/avgNum
                oNX.snr[oID] = SNR


        oNX.save()

def updateSNRSS(oFN):

        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['avgNumSS', 'numSignificantSequences', 'snrSS'])

        zeroCount = 0 #keeps tracks of the amount of simulations with 0 targets 
        for oID in oNX.snrSS:
                actualNum = oNX.numSignificantSequences[oID]
                avgNum = oNX.avgNumSS[oID]
         
                #give them one target for sims if none 
                if avgNum == 0:
                        zeroCount += 1
                        avgNum = 1
                
                SNR = float(actualNum)/avgNum
                oNX.snrSS[oID] = SNR

        print '# oRNA with 0 simulation targets: %s' % zeroCount 
        oNX.save()

if __name__ == "__main__":
        import sys
        #bioLibCG.submitArgs(updateAvgNumTargets, sys.argv)        
        #bioLibCG.submitArgs(updateSNR, sys.argv)        
        bioLibCG.submitArgs(updateAvgNumSS, sys.argv)        
        bioLibCG.submitArgs(updateSNRSS, sys.argv)        
