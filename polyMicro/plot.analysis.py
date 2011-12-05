import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast
import matplotlib.pyplot as plt

def plotBar(analysisFN):

    polyTable = cgNexusFlat.quickTable(('microSeq', 'string', '.', 0),
                                    ('readSeq', 'string', '.', 1),
                                    ('trimAmount', 'int', 0, 2),
                                    ('tailLength', 'int', 0, 3),
                                    ('letter', 'string', 'N', 4)) 
    polyNX = cgNexusFlat.Nexus(analysisFN, polyTable, ids = False)
    polyNX.load(['trimAmount', 'tailLength', 'letter'])

    letter_trim_tLen_count = {}

    for id in polyNX.ids:
        let = polyNX.letter[id]
        trim = polyNX.trimAmount[id]
        tLen = polyNX.tailLength[id]

        prevCount = letter_trim_tLen_count.setdefault(let, {}).setdefault(trim, {}).setdefault(tLen, 0)
        letter_trim_tLen_count[let][trim][tLen] = prevCount + 1


    letters = ["A", "T", "C", "G"]

    ind = range(1,5) 
    width = 0.35       # the width of the bars

    for trim in range(0, 8):
        for tLen in range(0,7):
            plotNum = (trim * 7) + tLen
            plt.subplot(8,7, plotNum)
            letAmount = []
            for let in letters: 
                letAmount.append(letter_trim_tLen_count[let].setdefault(trim, {}).setdefault(tLen, 0))
            #normalize
            theMax = max(letAmount)
            if theMax == 0:
                letAmount = [0,0,0,0]
            else:
                letAmount = [float(x+ .01)/theMax for x in letAmount]
            plt.bar(ind, letAmount, width)
            plt.title('%s,%s %s' % (trim, tLen, theMax) )
            plt.xticks([], [''])
            plt.yticks([], [''])

    plt.show()


    

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
