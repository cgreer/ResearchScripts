import bioLibCG
from cgNexus import Nexus
import cgDL
from cgAutoCast import autocast
from cgAutoKeyWord import autokey
from bioLibJA import subit

def checkDiFrequency(dataSeqsFN, simSeqsFN, outFN):

    dataSeqs = cgDL.listFromColumns(dataSeqsFN, [0], ['string'])
    simSeqs = cgDL.listFromColumns(simSeqsFN, [0], ['string'])
 
    def returnDiFreq(seqs):
        #collect di count
        di_freq = {}
        totalFrames = 0.0
        for seq in seqs:
            dis = bioLibCG.returnFrames(seq, 2)
            totalFrames += len(dis)
            for di in dis:
                di_freq[di] = di_freq.get(di, 0) + 1.0

        #convert to frequencey
        di_freq = dict( (x, di_freq[x]/totalFrames) for x in di_freq)
        return di_freq

    data_freqs = returnDiFreq(dataSeqs)
    sim_freqs = returnDiFreq(simSeqs)

    with open(outFN, 'w') as f:
        for di in data_freqs:
            f.write('%s\t%s\t%s\n' % (di, data_freqs.get(di, 0.0), sim_freqs.get(di, 0.0)))

        
    

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

