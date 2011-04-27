import bioLibCG
import cgDB
import cgOriginRNA
import matplotlib.pyplot as plt
from math import sqrt

def stdv(x): 
        n, mean, std = len(x), 0, 0 
        for a in x: 
                mean = mean + a 
        mean = mean / float(n) 
        for a in x: 
                std = std + (a - mean)**2
        std = sqrt(std / float(n-1))
        return std 

def totalSNR(oDir, filterList):
        
        fList = []
        f = open(filterList, 'r')
        for line in f:
                ls = line.strip().split('\t')
                fList.append(int(line.strip()))

        oDC = cgDB.dataController(oDir, cgOriginRNA.OriginRNA)
        id_oRNA = oDC.load()

        totalRun = 0 
        totalSim = 0
        simsTotals = []
        n = 0
        for oRNA in id_oRNA.values():
                if oRNA.id in fList:
                        if len(oRNA.filteredTargets) == 0:
                                continue
                        totalRun += len(oRNA.filteredTargets)
                        simsTotal = oRNA.avgNumSimulationTargets * 10
                        simsTotals.append(simsTotal)
                        totalSim += oRNA.avgNumSimulationTargets
                        n +=1
        
        print 'Total Number Targets for my run:', totalRun
        print 'Total Number Targets for Simulations:',  totalSim
        print 'SNR', float(totalRun)/float(totalSim)
        print 'stderr(%s)' % n, stdv(simsTotals)/sqrt(10)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(totalSNR, sys.argv)
