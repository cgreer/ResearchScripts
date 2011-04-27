import bioLibCG

def parseTargets(fN):

        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                micro = ls[3].split(',')
                micro.extend(ls[4].split(','))
                done = []
                for m in micro:
                        if m != 'None' and (m not in done):
                                print m
                                done.append(m)



if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(parseTargets, sys.argv)

