import bioLibCG
import cgNexusFlat

def cleanSet(fN, outFN):

        outF = open(outFN, 'w')
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                
                type1 = ls[13]

                if type1 == "processed_transcript_noncoding":
                        ls[15] = 'noncoding_gene'
                        outF.write('\t'.join(ls) + '\n')
                else:
                        outF.write(line)
                                        

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
