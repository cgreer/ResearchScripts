import bioLibCG

def process(fN):

        done = []
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                c1 = ls[0]
                if c1 not in done:
                        done.append(c1)
                else:
                        continue
                c = ls[0].split(':')
                gName = ls[1]

                print '%s\t%s\t%s\t%s' % (c[0], c[2], c[1], gName)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(process, sys.argv)
                
