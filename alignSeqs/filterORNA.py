import bioLibCG
import cgDB
import cgOriginRNA

def filterORNA(oDir, maxEndContig, maxTotalContig, minSNR, minNumTargets, keepDuplicates = False):
        if keepDuplicates == 'True': keepDuplicates = True
        if keepDuplicates == 'False': keepDuplicates = False

        maxEndContig, maxTotalContig = int(maxEndContig), int(maxTotalContig)
        minNumTargets = int(minNumTargets)
        minSNR = float(minSNR)

        oDC = cgDB.dataController(oDir, cgOriginRNA.OriginRNA)
        id_oRNA = oDC.load()

        for oRNA in id_oRNA.values():

                oRNA.passedFilter = True
                
                if len(oRNA.filteredTargets) < minNumTargets:
                        oRNA.passedFilter = False
                        cgOriginRNA.prettyPrint(oRNA, 'numTargets')
                        continue
                
                if oRNA.endContigLength > maxEndContig:
                        cgOriginRNA.prettyPrint(oRNA, 'endContig')
                        oRNA.passedFilter = False
                        continue

                if oRNA.totalContigLength > maxTotalContig:
                        cgOriginRNA.prettyPrint(oRNA, 'totalContig')
                        oRNA.passedFilter = False
                        continue

                if oRNA.snr < minSNR:
                        cgOriginRNA.prettyPrint(oRNA, 'SNR fail')
                        oRNA.passedFilter = False
                        continue

                if not keepDuplicates:
                        if oRNA.sequenceDuplicate:
                                cgOriginRNA.prettyPrint(oRNA, 'Duplicate Fail')
                                oRNA.passedFilter = False
                                continue
                
                print 'PASSED:', oRNA.id, ','.join(str(x) for x in oRNA.filteredTargets), oRNA.entropy, oRNA.avgNumSimulationTargets, oRNA.snr, oRNA.endContigLength, oRNA.sequence

        oDC.commit(id_oRNA)
        

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(filterORNA, sys.argv)
