import bioLibCG

def formatFile(fN):

        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split(',')
                transcript = ls[0][1:-1]
                miName = ls[1][1:-1]
                try:
                        location = int(ls[2])
                except:
                        continue

                if location == '':
                        continue

                print '%s\t%s\t%s' % (transcript, miName, location)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(formatFile, sys.argv)
                
