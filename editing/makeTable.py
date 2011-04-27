import bioLibCG
import cgEdit

def makeTable(fN, eFN):

        eSites = cgEdit.loadEditingSites(eFN)

        eID_eSite = {}
        for eSite in eSites:
                eID_eSite[eSite.ID] = eSite


        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                
                eID = int(ls[0])
                eSite = eID_eSite[eID]
                gName = ls[1]
                b = ls[8]
                a = ls[9]
                eRatio = eSite.eRatio
                eLoc = eSite.tcc

                print '%s\t%s\t%s\t%s\t%s' % (gName, eLoc, eRatio, b, a)


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(makeTable, sys.argv)
                
                        
