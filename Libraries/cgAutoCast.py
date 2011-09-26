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

def estimateType(var):
        '''guesses the str representation of the variable's type'''
        var = str(var) #important if the parameters aren't strings...
        for caster in (boolify, int, float, noneify):
                try:
                        return caster(var)
                except ValueError:
                        pass
        return var

def autocast(dFxn):
        def wrapped(*c, **d):
                cp = [estimateType(x) for x in c]
                dp = dict( (i, estimateType(j)) for (i,j) in d.items())
                return dFxn(*cp, **dp)
        
        return wrapped

def autoSig(dFxn):
        '''Does not work yet... just want to remind myself to do it.'''
        def wrapped(*c, **d):
                cp = [estimateType(x) for x in c]
                dp = dict( (i, estimateType(j)) for (i,j) in d.items())
                return dFxn(*cp, **dp)
        
        return wrapped


@autocast
def fxn2(one, two):
       print one, two, type(one), type(two)
        
if __name__ == "__main__":
        import sys

        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
