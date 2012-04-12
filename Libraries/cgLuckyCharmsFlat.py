'''Serialization Database'''
''' Do not have to worry about blanks with stringToX because db will fill in defaults for that'''

def stringToBool(s):
	'''Shorter words will decrease database size'''
	if s == 'T':
		return True
	else:
		return False

def boolToString(b):
	'''Shorter words will decrease database size'''
	if b:
		return 'T'
	else:
		return 'F'

def intListToString(l):
        if l:
                return ','.join([str(x) for x in l])
        else:
                return '.'

def stringToIntList(s):
	return [int(x) for x in s.split(',')]

def stringListToString(l):
        if l:
                return ','.join(l)
        else:
                return '.'

def stringToStringList(s):
        return s.split(',')

def boolListToString(l):
        if l:
                return ','.join([boolToString(x) for x in l])
        else:
                return '.'

def stringToBoolList(s):
	return [stringToBool(x) for x in s.split(',')]

def floatListToString(l):
        if l:
                return ','.join([str(x) for x in l])
        else:
                return '.'

def stringToFloatList(s):
	return [float(x) for x in s.split(',')]

def intListSpaceToString(l):
    if l:
        return ' '.join([str(x) for x in l])
    else:
        return '.'

def stringToIntListSpace(s):
        return [int(x) for x in s.split(' ')]

def stringListSpaceToString(l):
    if l:
        return ' '.join([str(x) for x in l])
    else:
        return '.'

def stringToStringListSpace(s):
        return [str(x) for x in s.split(' ')]

def returnSelf(s):
	return s

def getCasteFunction(dataType, fromDB = True):
	
	if dataType == 'int':
		if fromDB:
			return int
		else:
			return str

	if dataType == 'string':
		if fromDB:
			return returnSelf
		else:
			return returnSelf

	if dataType == 'bool':
		if fromDB:
			return stringToBool
		else:
			return boolToString
	
	if dataType == 'float':
		if fromDB:
			return float
		else:
			return str
	
	if dataType == 'intList':
		if fromDB:
			return stringToIntList
		else:
			return intListToString
	
	if dataType == 'stringList':
		if fromDB:
			return stringToStringList
		else:
			return stringListToString
		
	if dataType == 'boolList':
		if fromDB:
			return stringToBoolList
		else:
			return boolListToString
	
	if dataType == 'floatList':
		if fromDB:
			return stringToFloatList
		else:
			return floatListToString
	
        if dataType == 'stringListSpace':
		if fromDB:
			return stringToStringListSpace
		else:
			return stringListSpaceToString
	
        if dataType == 'intListSpace':
		if fromDB:
			return stringToIntListSpace 
		else:
			return intListSpaceToString 


