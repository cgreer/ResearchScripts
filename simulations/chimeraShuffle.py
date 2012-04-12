import bioLibCG
from cgNexus import Nexus
import cgDL
from cgAutoCast import autocast
from cgAutoKeyWord import autokey
from bioLibJA import subit
from cgLog import Logger 
logger = Logger(0)

import random

@autocast
def generateLeftRightChimera(oRNAFN, outFNBase):
    '''make map of 4mer --> left, 4mer rights
    NOTE: left and right do not include middle!''' 
    
    #FN is id/sequence
    sequences = cgDL.listFromColumns(oRNAFN, [1], ['string'])
     
    #update left/rights
    fourMer_left = {}
    fourMer_right = {}
    for seq in sequences:
        middle = seq[9:13]
        left, right = seq[:9], seq[13:]
        fourMer_left.setdefault(middle, set()).add(left)
        fourMer_right.setdefault(middle, set()).add(right)

    #lefts and rights mapped by middle 4
    for side in ['left', 'right']:
        outDict = eval('fourMer_%s' % side) #uh oh! EVAL TIME!
        with open(outFNBase + '.' + side, 'w') as f:
            for mer, seqs in outDict.iteritems():
                print mer, seqs
                f.write('%s\t%s\n' % (mer, ','.join(seqs)))

@autocast
def shuffleSingleSimulation(simFN, chimeraBaseFN, outFN, left = True): 
    '''left refers to which side you are KEEPING!  so left = True will only shuffle the right
    side'''

    #load chimera dictionaries
    middle_lefts = cgDL.dictFromColumns(chimeraBaseFN + '.left', (0, 1), ('string', 'stringList'))
    middle_rights = cgDL.dictFromColumns(chimeraBaseFN + '.right', (0, 1), ('string', 'stringList'))

    #given a sequence, make a chimera (make sure you're not using your own left/right)
    def shuffleOne(oSeq):
        oMiddle = oSeq[9:13]
        oLeft = oSeq[:9]
        oRight = oSeq[13:]

        lefts = middle_lefts[oMiddle][:]
        rights = middle_rights[oMiddle][:]

        lefts.remove(oLeft)
        rights.remove(oRight)

        if left:
            if not rights:
                return None
            shuffSeq = oSeq[:13] + random.choice(rights) 
        else:
            if not lefts:
                return None
            shuffSeq = random.choice(lefts) + oSeq[9:]

        return shuffSeq

    #shuffle each sequence in file
    with open(simFN, 'r') as f:
        with open(outFN, 'w') as fOut:
            with open(outFN + '.easyRead', 'w') as fOutER:
                for line in f:
                    id, oSeq = line.strip().split('\t')
                    chimeraSeq = shuffleOne(oSeq)
                    if not chimeraSeq: continue # NOTE: simulations with no shuffle will not be in output file
                    fOut.write('>%s\n%s\n\n' % (id, chimeraSeq))
                    fOutER.write('%s\t%s\n' % (id, chimeraSeq))
                        
if __name__ == "__main__":
    import sys
    assert sys.argv[1] in globals(), "Need name of fxn to run from command line!"
    fxnToRun = globals()[sys.argv[1]] 
    fxnToRun(*sys.argv[2:])

