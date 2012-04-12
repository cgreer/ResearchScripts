import bioLibCG
from cgNexus import Nexus
import cgDL
from cgAutoCast import autocast
from cgAutoKeyWord import autokey
from bioLibJA import subit
from cgLog import Logger 
logger = Logger(0)

from diShuffle import dinuclShuffle as diShuff
import re


def maskRunPositions(mask):
    '''get the positions of the runs of masked characters'''
    maskString = ''.join(['X' if x == True else 'O' for x in mask])
    runPositions = [(match.start(), match.end()) for match in re.finditer('X*', maskString)]
    runPositions = [(x,y) for x,y in runPositions if x!=y] #filter out 1-runs inside 2+ runs
    return runPositions

def shuffleMaskedSeq(seq, mask):
    
    #get masked run positions and sequences 
    runPositions = maskRunPositions(mask)
    runs = [seq[x:y] for x,y in runPositions] #y is already +1
    logger.log([runs, runPositions])
 
    #get the two letters that each run splits ("0" and "1" if at ends)
    def XCorrection(seq, position):
        return (seq[position - 1] if (position != 0) else "0")
    
    def YCorrection(seq, position):
        return (seq[position] if (position != len(seq)) else "1")#just len(seq) cuz y is +1 already

    runSplits = [(XCorrection(seq, x), YCorrection(seq, y)) for x,y in runPositions]
    logger.log(runSplits)

    #get sequence without masked parts (to be shuffled)
    shuffleSeq = seq[:]
    for run in runs:
        shuffleSeq = shuffleSeq.replace(run, '')


    #shuffle seq
    shuffleSeq = diShuff(shuffleSeq)

    #put masked sequences back into shuffled seqs (diseqs should still be there)
    finalShuffledSeq = shuffleSeq
    for run, runSplit in zip(runs, runSplits):
        if runSplit[0] == "0":
            finalShuffledSeq = run + finalShuffledSeq
        elif runSplit[1] == "1":
            finalShuffledSeq += run
        else:
            splitDi = '%s%s' % runSplit
            splitPosition = finalShuffledSeq.find(splitDi)
            finalShuffledSeq = finalShuffledSeq[:splitPosition + 1] + run + finalShuffledSeq[splitPosition + 1:]

    return finalShuffledSeq


if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

