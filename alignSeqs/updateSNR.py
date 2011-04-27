import bioLibCG
import cgDB
import cgOriginRNA

def updateAvgNumTargets(oDir):
       
        oID_numTargets = {}

        for i in range(0,10):
                print i
                simDirRNA = '/home/chrisgre/scripts/simulations/simsk50Filtered/simulation.%s/oRNA' % i
                oDC = cgDB.dataController(simDirRNA, cgOriginRNA.OriginRNA)
                id_sRNA = oDC.load()

                for id, sRNA in id_sRNA.items():
                        currTargets = oID_numTargets.get(id, 0)
                        oID_numTargets[id] = currTargets + len(sRNA.filteredTargets)
        
        #now save it
        oDC = cgDB.dataController(oDir, cgOriginRNA.OriginRNA)
        id_oRNA = oDC.load()

        for oRNA in id_oRNA.values():
                totalNum = oID_numTargets.get(oRNA.id, 0)
                avgNum = float(totalNum)/float(10.0)
                oRNA.avgNumSimulationTargets = avgNum

        oDC.commit(id_oRNA)

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
                oDC = cgDB.dataController(simDirRNA, cgOriginRNA.OriginRNA)
                id_sRNA = oDC.load()
                
                groupTotal = 0
                for id, sRNA in id_sRNA.items():
                        if not id in fList:
                                continue
                        groupTotal +=  len(sRNA.filteredTargets)
                print groupTotal                        
        

def updateSNR(oDir):

        oDC = cgDB.dataController(oDir, cgOriginRNA.OriginRNA)
        id_oRNA = oDC.load()

        for oRNA in id_oRNA.values():
                actualNum = len(oRNA.filteredTargets)
                avgNum = oRNA.avgNumSimulationTargets
                
                if avgNum == 0: avgNum = .01
                
                SNR = float(actualNum)/avgNum
                oRNA.snr = SNR


        oDC.commit(id_oRNA)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(updateAvgNumTargets, sys.argv)        
        bioLibCG.submitArgs(updateSNR, sys.argv)        
        #bioLibCG.submitArgs(countForError, sys.argv)        
