import bioLibCG
import cgGenes3
import dumpObj


def doit(fN):
        geneSet = cgGenes3.createGeneSetEditing(fN)

        for transcript in geneSet.transcripts:
                if transcript.id == 'NM_031422':
                        dumpObj.dumpObj(transcript)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(doit, sys.argv)
                        



