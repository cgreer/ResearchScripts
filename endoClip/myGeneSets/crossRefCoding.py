import bioLibCG
import cgNexusFlat


def crossRefDict(fN, keyCol, valCol, keyCasteFxn, valCasteFxn):
        
        key_val = {}
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                k, v = keyCasteFxn(ls[keyCol]), valCasteFxn(ls[valCol])
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
                        ls.append('NA')

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

def crossRefCoding(tFN, infoFN):

        #make cref dict
        tID_cStatus = crossRefDict(infoFN, 0, 1, str, str)
        crossUpdate(tFN, 0, tID_cStatus)

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
