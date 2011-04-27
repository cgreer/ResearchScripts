import bioLibCG
import cgEdit

def getEditInfo(fN, idList):
        
        eSites = cgEdit.loadEditingSites(fN)
        
        idDict = {}
        for eSite in eSites:
                idDict[eSite.ID] = eSite

        list = []
        f = open(idList, 'r')
        for line in f:
                ls = line.strip().split('\t')
                list.append(int(ls[0]))

        for id in list:
                eSite = idDict[id]
                print eSite.ID, '%s:%s' % (eSite.chromosome, eSite.coordinate), eSite.gene, eSite.eRatio
               
if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(getEditInfo, sys.argv)
               
