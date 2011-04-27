import bioLibCG
import cgDB
import cgOriginRNA
import matplotlib.pyplot as plt
import math

def signalHistoOLD(resultsFN, avgTargetFN):

        id_avgNum = {}
        f = open(avgTargetFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = int(ls[0])
                avgNum = float(ls[1])
                id_avgNum[id] = avgNum

        histVals = []
        f = open(resultsFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = int(ls[0])
                numTargets = len(ls[4].split(','))
                try:
                        avgNum = id_avgNum[id]
                except KeyError:
                        avgNum = .01
                sigNoiseRatio = numTargets/avgNum
                histVals.append(math.log(sigNoiseRatio, 2))
                #histVals.append(sigNoiseRatio)

        plt.hist(histVals, 30, facecolor='b', alpha = .75)
        plt.xlabel('Signal/Noise')
        plt.ylabel('Number of Origin RNAs')
        plt.show()

def signalHisto(oDir, title = 'SNR'):

        oDC = cgDB.dataController(oDir, cgOriginRNA.OriginRNA)
        id_oRNA = oDC.load()

        histVals = []
        for oRNA in id_oRNA.values():
                if not oRNA.passedFilter:
                        continue

                histVals.append(math.log(oRNA.snr, 2))

        plt.hist(histVals, 30, facecolor='b', alpha = .75)
        plt.axis([-4,10,0,30])        
        plt.title('%s' % title)
        plt.xlabel('log2(Signal/Noise)')
        plt.ylabel('Number of Origin RNAs')
        plt.show()
                        


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(signalHisto, sys.argv)
