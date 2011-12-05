import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast


#the var must be passed to lambda because it becomes a local variable inside lamda.  Can put local variables into global namespace?

def makeFxn(expr, g = globals(), l = locals()):
   
    print 'globals', g
    print 'loc', l
    return eval('lambda a=a: %s' % expr, g, l)

def parseExpression(expr):
    #this has two variables, a function, and a literal
    $entropy > fxn(a) + 1 + 'aaa'



def callFxn(expr):
   
    a = 3
    fxn = makeFxn(expr, globals(), locals())

    return fxn()

def returnOne():
    return 1

print callFxn('returnOne() + 3 + a')
