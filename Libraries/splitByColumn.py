import bioLibCG

def splitByColumn(fN, column, splitModifier = None, outDir = ""):
        column = int(column)
        if splitModifier == "None":
            splitModifier = None

        cValue_file = {}

        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                
                if splitModifier == 'chromStrand':
                        columnValue = '%s.%s' % (ls[column].split(':')[0], ls[column].split(':')[1])
                elif splitModifier == 'editing':
                        columnValue = '%s.%s' % (ls[1], ls[15])
                elif splitModifier == 'geneDefinition':
                        columnValue = '%s.%s' % (ls[1], ls[2])
                elif splitModifier == 'oneNT':
                        columnValue = '%s' % (ls[column].split(':')[0])
                elif splitModifier == 'geneSpan':
                        columnValue = '%s.%s' % (ls[column], ls[column + 1])
                else:
                        columnValue = ls[column]
                
                if columnValue not in cValue_file:
                    if outDir:
                        cValue_file[columnValue] = open(outDir + '/' + fN.strip().split('/')[-1] + '.' + columnValue, 'w')
                    else:
                        cValue_file[columnValue] = open(fN + '.' + columnValue, 'w')

                cValue_file[columnValue].write(line)

        #close all files
        for cValue, file in cValue_file.iteritems():
                file.close()


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(splitByColumn, sys.argv)
