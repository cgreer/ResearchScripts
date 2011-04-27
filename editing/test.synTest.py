import bioLibCG
import cgGenes3
import cgSeqMod

def test(gFN):
        
        map = cgSeqMod.loadCodonMap('hg19')
        geneSet = cgGenes3.createGeneSetEditing(gFN)

        for transcript in geneSet.transcripts:
                if '_coding' in transcript.tType:
                        try:
                                print transcript.id
                                print transcript.getMRNA(coding = True)
                                print cgSeqMod.translateRNA(transcript.getMRNA(coding = True), map)
                        except:
                                print 'fail'

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(test, sys.argv)
