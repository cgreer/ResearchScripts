import GenomeFetch as gf
import bioLibCG

def quickSeq(tcc, assembly):

        myG = gf.GenomeFetch(assembly)
        print myG.getSequence(tcc)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(quickSeq, sys.argv)

