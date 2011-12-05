import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast
from miRTable import miR 
from cgAlignmentFlat import cgAlignment

def revMap(NX, attName):
    '''get reverse mapping, always key --> list'''

    att = eval('NX.attName')

    att_ids = {}
    for id in NX.ids:
        att_ids.setdefault(att[id], []).append(id)

    return att_ids        
    

def updateClippedSeq(mFN):
    
    NX = cgNexusFlat.Nexus(mFN, miR)
    NX.load(['sequence', 'clippedSeq'])

    for id in NX.ids:

        NX.clippedSeq[id] = NX.sequence[id][:-4]

    NX.save()        

def tabIt(*args):

    return '\t'.join([str(x) for x in args])

def updatePolySeqs(mFN, readsFN, alignFN):

    tim = bioLibCG.cgTimer()
    tim.start()
    variousAs = ["A" * x for x in range(1,20)]
    variousGs = ["G" * x for x in range(1,20)]
    variousTs = ["T" * x for x in range(1,20)]
    variousCs = ["C" * x for x in range(1,20)]

    letter_variousLetters = [ ("A", variousAs),
                            ("G", variousGs),
                            ("T", variousTs),
                            ("C", variousCs)]


    checkRange = range(1,8)

    NX = cgNexusFlat.Nexus(mFN, miR)
    NX.load(['sequence', 'polySeqs'])
    #print 'load micro', tim.split() 

    reads = cgNexusFlat.quickTable(('read','string', '.', 1))
    rNX = cgNexusFlat.Nexus(readsFN, reads)
    rNX.load(['read'])
    #print 'load reads', tim.split() 

    aNX = cgNexusFlat.Nexus(alignFN, cgAlignment)
    aNX.load(['sID', 'tID'])
    #print 'load alignments', tim.split() 

    for id in aNX.ids:

        theRead = rNX.read[aNX.sID[id]]
        mID = aNX.tID[id]
        microSeq = NX.sequence[mID]

        #may be a read for expression, but wont count...
        if theRead in microSeq: continue

        #just for expression
        if microSeq == theRead: 
            print tabIt(microSeq, theRead, 0, 0, "N")

        #first check full
        elif microSeq in theRead and (len(theRead) != len(microSeq)):
            tail = theRead.split(microSeq)[1]
            for let, variousLetters in letter_variousLetters:
                if tail in variousLetters:
                    print tabIt(microSeq, theRead, 0, len(tail), let)

        #now check trimmed (cant do [:-0])
        else:
            for i in checkRange:
                if microSeq[:-i] in theRead and (len(theRead) != len(microSeq[:-i])):
                    tail = theRead.split(microSeq[:-i])[1]
                    for let, variousLetters in letter_variousLetters:
                        if tail in variousLetters:
                            print tabIt(microSeq, theRead, i, len(tail), let)
                            print "TRIMMED"
                    break #dont trim after the first trimmed one works                           




if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
