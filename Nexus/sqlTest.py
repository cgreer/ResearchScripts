import sqlAddOn
import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast
import cgOriginRNAFlat

def testSQL(oFN):
    
    NX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
    NX.load(['entropy', 'filteredTargets'])


    sNX = sqlAddOn.select(NX, select = ['filteredTargets', 'entropy'],
                            where = '$entropy > 1.7')

    
    ents = sqlAddOn.collect(sNX, select = ['filteredTargets'],
                                where = 'len( $filteredTargets ) > 1',
                                unique = False)
    global a
    a = 3
    sqlAddOn.update(sNX, sets = ['$entropy = $entropy + a',
                                '$filteredTargets = $filteredTargets'], g = globals())
    print sNX

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
