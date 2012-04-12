import bioLibCG
from cgNexus import Nexus
import cgDL
from cgAutoCast import autocast
from cgAutoKeyWord import autokey
from bioLibJA import subit
from cgLog import Logger 
logger = Logger(0)

def singleView(oFN, unCleanOFN, aFN, dFN):

    oID_aID_dID = {}
    oID_siblingIDs = {}
    with open(oFN, 'r') as f:
        for line in f:
            ls = line.strip().split('\t')
            oID = int(ls[0])
            oID_siblingIDs[oID] = [int(x) for x in ls[19].split(',')]

            aIDs = [int(x) for x in ls[8].split(',')]
            oID_aID_dID[oID] = dict( (aID, None) for aID in aIDs)

    #lines have to come from original unclean file to get 
    #simulation lines as well
    oID_line = {}
    with open(unCleanOFN, 'r') as f:
        for line in f:
            ls = line.strip().split('\t')
            oID = int(ls[0])
            oID_line[oID] = line

    aID_line = {}
    dIDs = set()
    with open(aFN, 'r') as f:
        for line in f:
            ls = line.strip().split('\t')
            aID, oID, dID = [int(x) for x in ls[0:3]]
            if aID in oID_aID_dID.get(oID, []):
                aID_line[aID] =  " " + line
                oID_aID_dID[oID][aID] = dID
                dIDs.add(dID)
    
    dID_line = {}
    with open(dFN, 'r') as f:
        for line in f:
            ls = line.strip().split('\t')
            dID = int(ls[0])
            if dID in dIDs: 
                dID_line[dID] =  "  " + line

    #now write them into a single view
    for oID in oID_aID_dID:
        print 
        #oRNAs
        print oID_line[oID],
        for sID in oID_siblingIDs[oID]:
            if sID == oID: continue
            print oID_line[sID],
        for aID in oID_aID_dID[oID]:
            dID = oID_aID_dID[oID][aID]
            print aID_line[aID],
            print dID_line[dID],

if __name__ == "__main__":
    import sys
    assert sys.argv[1] in globals(), "Need name of fxn to run from command line!"
    fxnToRun = globals()[sys.argv[1]] 
    fxnToRun(*sys.argv[2:])

