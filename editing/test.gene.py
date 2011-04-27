import bioLibCG
import cgGenes3
import cgSeqMod

def testit(gFN):

        geneSet = cgGenes3.createGeneSetEditing(gFN)

        map = cgSeqMod.loadCodonMap('hg19')
        for gene in geneSet.genes:
                for transcript in gene.transcripts:
                        try:
                                print ''
                                mRNA = transcript.getMRNA(coding = True)
                                i = transcript.getRelativePositionMRNA(35872409)
                                if i == -1:
                                        continue
                                print transcript.id
                                print i
                                print mRNA[:i], mRNA[i], mRNA[i + 1:]
                                print cgSeqMod.translateRNA(mRNA, map)
                        except:
                                pass


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(testit, sys.argv)

