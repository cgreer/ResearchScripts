import bioLibCG
from cgNexusFlat import Field
                        

class MicroSeq:

        seq = Field('string', '.', 1)
        longestMono = Field('int', 0, 2)
        longestDi = Field('int', 0, 3)
        longestTri = Field('int', 0, 4)
        longestQuad = Field('int', 0, 5)
    

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
