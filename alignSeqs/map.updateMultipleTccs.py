import cgDB
import cgOriginRNA
import bioLibCG

def updateMultipleTccs(oDir, mappedFN):

        #parse bowtie
        f = open(mappedFN, 'r')
        oID_tccs = {}
        for line in f:
                ls = line.strip().split('\t')
                oID = int(ls[0])
                strand, chrom, firstCoord = bioLibCG.switchStrandFormat(ls[1]), ls[2], int(ls[3])
                secondCoord = firstCoord + len(ls[4]) - 1
                tcc = bioLibCG.makeTcc(chrom, strand, firstCoord, secondCoord)
                
                oID_tccs.setdefault(oID, []).append(tcc)


        oDC = cgDB.dataController(oDir, cgOriginRNA.OriginRNA)
        id_oRNA = oDC.load()

        for id, oRNA in id_oRNA.items():

                oRNA.tccs = oID_tccs[id]
        
        oDC.commit(id_oRNA)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(updateMultipleTccs, sys.argv)
