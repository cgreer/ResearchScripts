import bioLibCG

def boolify(s):

    if s == 'True' or s == 'true':
        return True
    if s == 'False' or s == 'false':
        return False
    raise ValueError('Not Boolean Value!')

def noneify(s):
    ''' for None type'''
    if s == 'None': return None
    raise ValueError('Not None Value!')

def listify(s):
    '''will convert a string representation of a list
    into list of homogenous basic types.  type of elements in
    list is determined via first element and successive 
    elements are casted to that type'''

    #this cover everything?
    if "," not in s:
        raise ValueError('Not a List')

    #derive the type of the variable
    loStrings = s.split(',')
    elementCaster = None
    for caster in (boolify, int, float, noneify, str):
        try:
            caster(loStrings[0])
            elementCaster = caster
            break
        except ValueError:
            pass

    #cast all elements
    try:
        castedList = [elementCaster(x) for x in loStrings]
    except ValueError:
        raise TypeError("Autocasted list must be all same type")

    return castedList


def estimateTypedValue(var):
    '''guesses the str representation of the variable's type'''

    #dont need to guess type if it is already un-str typed (not coming from CLI)
    if type(var) != type('aString'):
        return var 

    #guess string representation, will default to string if others dont pass
    for caster in (boolify, int, float, noneify, listify, str):
        try:
            return caster(var)
        except ValueError:
            pass

def autocast(dFxn):
    def wrapped(*c, **d):
        cp = [estimateTypedValue(x) for x in c]
        dp = dict( (i, estimateTypedValue(j)) for (i,j) in d.items())
        return dFxn(*cp, **dp)
    return wrapped

@autocast
def test(one, two):
       print one, two, type(one), type(two)
        
if __name__ == "__main__":
    import sys
    assert sys.argv[1] in globals(), "Need name of fxn to run from command line!"
    fxnToRun = globals()[sys.argv[1]] 
    fxnToRun(*sys.argv[2:])
