import bioLibCG
import cgNexusFlat
import cgDL
import inspect
from cgAutoCast import autocast

def continuityCheck(theList):

    for i,item in enumerate(theList):
        if i == 0 or i == len(theList) - 1:
            continue

        if theList[i - 1] + 1 == item:
            return True

        if theList[i + 1] - 1 == item:
            return True

    return False

def autokey(dFxn):

    def wrapped(*nonkw, **kw):
        '''parameters from command line will be seen as non-kw (in c)
        check for duplicate key words
        check for keywords without following parameters'''
        
        #get possible keywords, positions of keywords
        fxnArgNames = inspect.getargspec(dFxn)[0]
        possibleKeywords = ['-%s' % x for x in fxnArgNames]
        keyPositions = [i for i, name in enumerate(nonkw) if name in possibleKeywords]

        #no keywords, just return
        if not keyPositions:
            return dFxn(*nonkw, **kw)

        #error checks
        if len(nonkw) > (keyPositions[-1] + 2):
            raise NameError("non-keyword arguments can not follow kw args")

        if continuityCheck(keyPositions):
            raise NameError("two keyword designations in a row!")

        for theKW in possibleKeywords:
            if nonkw.count(theKW) > 1:
                raise NameError("keyword was used twice!")

        #update nonkw
        newNonkw = [x for x in nonkw[:keyPositions[0]] ]

        #update kw
        for position in keyPositions:
            key, val = nonkw[position][1:], nonkw[position + 1] #1 on index for '-'
            kw[key] = val

        return dFxn(*newNonkw, **kw)
    
    return wrapped
    

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

