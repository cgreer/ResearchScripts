import bioLibCG

def splitByColumn(fN, column, splitModifier = None):
        column = int(column)

        cValue_file = {}

        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                
                if splitModifier == 'chromStrand':
                        columnValue = '%s.%s' % (ls[column].split(':')[0], ls[column].split(':')[1])
                else:
                        columnValue = ls[column]
                
                if columnValue not in cValue_file:
                        cValue_file[columnValue] = open(fN + '.' + columnValue, 'w')

                cValue_file[columnValue].write(line)

        #close all files
        for cValue, file in cValue_file.iteritems():
                file.close()


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(splitByColumn, sys.argv)
