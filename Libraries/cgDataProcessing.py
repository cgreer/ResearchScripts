import bioLibCG
        
def crossRefDict(fN, keyCol, valCol, keyCasteFxn, valCasteFxn):
        
        key_val = {}
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                try:
                        k, v = keyCasteFxn(ls[keyCol]), valCasteFxn(ls[valCol])
                except IndexError:                        
                        k, v = keyCasteFxn(ls[keyCol]), valCasteFxn(None)
                key_val[k] = v

        f.close()
        return key_val


def crossUpdate(uFN, primaryKey, crossDict):

        f = open(uFN, 'r')
        newLines = []
        for line in f:
                ls = line.strip().split('\t')
                primaryID = ls[primaryKey]
                try:
                        ls.append(crossDict[primaryID])
                except KeyError:
                        ls.append('.')

                newLines.append('\t'.join([str(x) for x in ls]) + '\n')
        f.close()

        f = open(uFN, 'w')
        f.writelines(newLines)
        f.close()

def cgJoin(uFN, jFN, updateKeyColumn, joinKeyColumn, joinValColumn):
        '''add parameter from one column of joined file to updated File
        the values compared are strings.'''
        updateKeyColumn, joinKeyColumn, joinValColumn = int(updateKeyColumn), int(joinKeyColumn), int(joinValColumn)
        key_parameter = crossRefDict(jFN, joinKeyColumn, joinValColumn, str, str)
        crossUpdate(uFN, updateKeyColumn, key_parameter)

def removeColumn(fN, columnNumber, blank = False):
        if blank == 'True':
                blank = True
        elif blank == 'False':
                blank = False
        columnNumber = int(columnNumber)                

        f = open(fN, 'r')
        newLines = []
        for line in f:
                ls = line.strip().split('\t')
                
                if blank:
                        ls[columnNumber] = '.'
                else:
                        del ls[columnNumber]

                newLines.append('\t'.join(ls) + '\n')
        f.close()

        f = open(fN, 'w')
        f.writelines(newLines)
        f.close()

def addColumn(fN, afterColumnNumber, addText):
        afterColumnNumber = int(afterColumnNumber)

        f = open(fN, 'r')
        newLines = []
        for line in f:
                ls = line.strip().split('\t')
                
                ls_1 = ls[:afterColumnNumber + 1]
                ls_2 = ls[afterColumnNumber + 1:]

                ls_1.append(addText)
                ls = ls_1 + ls_2
                newLines.append('\t'.join(ls) + '\n')
                
        f.close()

        f = open(fN, 'w')
        f.writelines(newLines)
        f.close()


if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

