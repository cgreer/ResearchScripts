import bioLibCG


def uniqueColumn(fN, column, whole = False):

        u = {}
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                u[ls[int(column)]] = line.strip()

        if whole:
                for i in u:
                        print u[i]
        else:
                for i in u:
                        print i

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(uniqueColumn, sys.argv)
                        

