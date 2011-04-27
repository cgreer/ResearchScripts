import bioLibCG
import cgDB
import cgOriginRNA

def plotContextPie(oDir, contextFN):
        
        oDC = cgDB.dataController(oDir, cgOriginRNA.OriginRNA)
        id_oRNA = oDC.load()

        oID_tTypes = {}
        f = open(contextFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = int(ls[0])
                tType = ls[5]
                oID_tTypes.setdefault(id, []).append(tType)


        tType_count = {}
        for oRNA in id_oRNA.values():
                if not oRNA.passedFilter:
                        continue
                

                for tType in oID_tTypes[oRNA.id]:
                        num = tType_count.get(tType, 0)
                        tType_count[tType] = num + 1

        for tType, count in tType_count.items():
                print tType, count



if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(plotContextPie, sys.argv)
