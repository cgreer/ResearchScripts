import bioLibCG
import cgOriginRNAFlat
import matplotlib.pyplot as plt
from math import sqrt
import cgNexusFlat

def stdv(x): 
        n, mean, std = len(x), 0, 0 
        for a in x: 
                mean = mean + a 
        mean = mean / float(n) 
        for a in x: 
                std = std + (a - mean)**2
        std = sqrt(std / float(n-1))
        return std 

def totalSNR(oFN):
        
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['snr', 'filteredTargets', 'avgNumSimulationTargets', 'passedFilter'])

        totalRun = 0 
        totalSim = 0
        simsTotals = []
        n = 0
        highSNR = 0
        for oID in oNX.avgNumSimulationTargets:
                
                #filter out
                if not oNX.passedFilter[oID]:
                        continue 
                           
                #collect stats
                totalRun += len(oNX.filteredTargets[oID])
                simsTotal = oNX.avgNumSimulationTargets[oID] * 10
                simsTotals.append(simsTotal)
                totalSim += oNX.avgNumSimulationTargets[oID]
                
                if oNX.snr[oID] > 2:
                        highSNR += 1
                n +=1
       
        print oFN
        print 'Total Number Targets for my run:', totalRun
        print 'Total Number Targets for Simulations:',  totalSim
        print 'SNR', float(totalRun)/float(totalSim)
        print 'Total oRNA:', n, 'Total oRNA w/ SNR > 2:', highSNR
        print '\n' 
        #print 'stderr(%s)' % n, stdv(simsTotals)/sqrt(10)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(totalSNR, sys.argv)
