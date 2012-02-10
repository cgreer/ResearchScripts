import bioLibCG
from cgAutoCast import autocast
import cgLuckyCharmsFlat

def dictFromColumns(fN, keyValCols, keyValTypes, assumeUnique = True):
    '''Cols and Types are 2-tuple with key first, val second
    assumeUnique will raise Error True and non-unique value arises'''

    keyCasteFxn = cgLuckyCharmsFlat.getCasteFunction(keyValTypes[0])
    valCasteFxn = cgLuckyCharmsFlat.getCasteFunction(keyValTypes[1])
    keyCol = keyValCols[0]
    valCol = keyValCols[1]

    key_val = {}
    f = open(fN, 'r')
    for line in f:
        ls = line.strip().split('\t')
       
        key, val = keyCasteFxn(ls[keyCol]), valCasteFxn(ls[valCol])

        if assumeUnique:
            if key in key_val:
                raise NameError("Mapping is not 1 to 1")
            else:
                key_val[key] = val

        else:
            key_val.setdefault(key, []).append(val)
    f.close()
   
    return key_val

def listFromColumns(fN, columns, valTypes, mergeType = 'lol'):
    '''multiple columns should either go into multiple lists (list of lists)
    or merge into the same list
    mergeType = lol or merge (list of list or merge)'''

    #checks
    if len(columns) != len(valTypes):
        raise NameError("Must provide Types for ALL columns")
  
    
    colNum_casteFxn = dict( (i, cgLuckyCharmsFlat.getCasteFunction(valTypes[i])) for i in range(len(columns)) )

    lol = [list() for i in range(len(columns))] 
    f = open(fN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        for i in columns:
            lol[i].append(ls[i])
    f.close()
    
    if len(lol) == 1:
        return lol[0]
    elif mergeType == 'lol':
        return lol
    elif mergeType == 'merge':
        mergedList = []
        [mergedList.extend(x) for x in lol]
        return mergedList
    else:
        raise NameError("WTH?!")
        

    
    
    
if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

