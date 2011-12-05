import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast

# can use & for "other" NX.  I think sql statements only have two tables at a time...no but mine will
# e.g., '$sequence == &sequence'
#may have to use locals
#everything you do is basically a left join (all left and part of right if present)
#you can even put custom fxns in expression I believe

#Note: passing more parameters to a fxn costs time...

def createConditionalFxn(condExpression, g = globals()):
    '''$att is for firstNX, &att is for 2nd.  Only two supported
    The dynamic variables ($/&)att replaced with fNX.att[x]
    Each condFxn has 4 variables even if only 2 are used
    there must be spaces in expression el can't parse.  This
    can be fixed by checking for att in class...'''
    
   
    for dynMark, NX, idVar in ( ('$', 'fNX.', '[x]'), ('&', 'sNX.', '[y]') ):
        #get dynamic features
        dynVars = [x for x in condExpression.split() if x.startswith(dynMark)]
        dynVars = list(set(dynVars))
        
        #replace dynVars with NX.att[$/&]
        for var in dynVars:
            condExpression = condExpression.replace(var, NX + var[1:] + idVar)
   
    print 'a' in g
    lambdaText = 'lambda fNX, x, sNX, y: %s' % condExpression
    print lambdaText
    return eval(lambdaText, g)
     
def collect(NX, select = [], where = '', unique = False):
    '''put items in column selected columns into a list
    if colType is list --> add all items
    if >1 column selected --> add items in both columns into one list'''
    

    #create where fxn
    whereFxn = createConditionalFxn(where)

    #occupy selections to collect
    listType = type([])
    selectAttributes = [] #(type, attribute) 
    for sel in select:
        att = getattr(NX, sel)
        theType = type(att[att.iterkeys().next()]) #first item in att (gen is much faster)
        selectAttributes.append( (theType, att ) )

    if unique:
        collection = set()
    else:        
        collection = []
    
    for id in NX.ids:

        #Does it pass where statement?
        if not whereFxn(NX, id, 0, 0): continue

        #do collection
        for theType, att in selectAttributes:
            if theType == listType:
                if unique:
                    for i in att[id]:
                        collection.add(i)
                else:
                    collection.extend(att[id])
            else:
                if unique:
                    collection.add(att[id])
                else:    
                    collection.append(att[id])
                     
    return collection         
                    
def select(NX, select = [], where = '', retain = False):
    '''select columns from NX and place into a new NX

    Saving: TODO 
    If you select a subset of lines and change some 
    values in the line, when you save it will still
    write every line in the original packet...

    Retaining table structure is optional (def false)
    because if you want to filter whole db then just
    select all attributes'''

    # "*" selection
    if select[0] == '*':
        select = NX._selectedAttNames

    #create empty NX to fill
    #TODO create a separate save property FN
    selectNX = cgNexusFlat.shellNexus(NX, NX._dataFileName + '.selected', select)
    selectNX.initializeMasterDict()
    selectNX.bindAttributes(select)
    selectNX.linkIDsToColumn()
    
    #collapse columns if table structure isn't retained
    if not retain:
        selectNX.collapseColumnNumbers(select)

    #create where fxn
    whereFxn = createConditionalFxn(where)

    #occupy selections to collect
    NXAttributes = []
    selectAttributes = [] #(type, attribute) 
    for sel in select:
        selectAttributes.append(getattr(selectNX, sel))
        NXAttributes.append(getattr(NX, sel))
    bothAttributes = zip(NXAttributes, selectAttributes)

    #fill selection
    for id in NX.ids:

        if not whereFxn(NX, id, 0, 0): continue

        #fill new NX
        for NXAtt, selectAtt in bothAttributes:
            selectAtt[id] = NXAtt[id]

    return selectNX            

def join(NX, NX2, select = [], where = '', on = '', into = []):
    ''' merge two NXs
    Duplicate Lines:
    Keys are unique so there might be some problems for many to many.
    I think you can invert the join direction to make many to one
    If there is more than one match in 2nd NX the column should be set 
    to list type.  The join should automatically add the list...
    This list functionality may take care of o-to-m and m-to-m

    Select is for 2nd NX attributes.  All attributes from the first NX
    will carry over (can do a selection if you want to narrow in a sep step).

    on defaults to '$id == &id' which would be the 0th columns of each file.
    Have to make custom fork for that.  Else on will make dict of 2nd NXs
    values as keys against the id of that line.  {primkey: lineID} or 
    {primKey: [lineIDs]} for list types.  Non one-to-one mappings MUST
    have list type specified or a NameMapError should occur on failed append.

    Explicitly mapped column names in into ('&seq', '$sequence')
    are used before default (same name in two tables) If explicit and
    default don't cover all names --> raise NameMapError
    
    The above rules mean that left joins cannot be made on NXs that don't have
    proper tables made already... perhaps make optional variable?  Allow dynamic
    alteration of table structure?
    '''

    pass

def update(NX, sets = [], where = '', g = {}):
    ''' set is '$sequence = 'AAAA' .  Might want to check
    for equals sign and only one attribute on left side
        
    Gonna need another condition fxn generator for Update
    Should put all sets in one lambda to save time...
    
    TODO: $prop = $prop.append(2) doen not work cuz the expr 
    ends up being fNX.prop.append(2)[x]...fix condExpression

    lambda fxns cannot set variables, only return expressions'''

    #validate sets
    for s in sets:
        if s.count('=') != 1 or s.split('=')[0].count('$') != 1:
            raise NameError("Update Set Is Improper")

    if where:
        whereFxn = createConditionalFxn(where)        


    #establish binding attributes and their update fxns
    attributes, expressions = [], []
    for set in sets:
        
        #parse strings for expression and attribute
        attName = set.split('=')[0].strip()[1:]
        expressionString = set.split('=')[1]
        print attName, expressionString

        attExpression = createConditionalFxn(expressionString, g = g)
        attribute = eval('NX.%s' % attName )
        attributes.append(attribute)
        expressions.append(attExpression)

    att__expression = zip(attributes, expressions)

    for id in NX.ids:

        #check whereFxn
        if where:
            if not whereFxn(NX, id, 0, 0): continue

        #update
        for att, expr in att__expression:
            att[id] = expr(NX, id, 0, 0)

