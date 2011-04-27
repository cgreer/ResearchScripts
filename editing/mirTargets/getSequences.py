import bioLibCG

def getSeqs(fN):

        name_seq = {}

        f = open(fN, 'r')
        mirName = None
        for line in f:

                
                if line.startswith('>'):
                        if not 'Homo sapiens' in line:
                                ls = line.strip().split(' ')
                                mirName = ls[-1]
                else:
                                
                        seq = line.strip()

                        name_seq[mirName] = seq




        f.close()

        for name in name_seq:
                print '%s\t%s' % (name, name_seq[name])



if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(getSeqs, sys.argv)
                

