import math

eVal = math.e

def getLam(distList):
	'''given a list of points from a distribution, find the lambda of a exponential
	distibution'''
	sum = math.fsum(distList)		
	mean = sum/float(len(distList))
	if mean == 0:
		return float('inf')
	else:
		return float(1)/mean

def getPValExp(x, lam):
	x = float(x)
	lam = float(lam)
	lamx = lam*x
	try:
		pVal = float(1)/(eVal**lamx)
	except OverflowError:
		print 'pVal too small'
		pVal = .00000001
		
	
	return pVal
	
