#Take a file and sort it by XX
from operator import itemgetter
numRange = range(1,99)
letRange = ['M', 'X', 'Y', 'Z']
import subprocess

def returnChromToNumDict():
	chromToNum = {}
	for i in numRange:
		chrom = 'chr%s' % i
		chromToNum[chrom] = i
	for i, let in enumerate(letRange):
		chrom = 'chr%s' % let
		chromToNum[chrom] = i + 1000
	return chromToNum

def returnNumToChromDict():
	numToChrom = {}
	for i in numRange:
		chrom = 'chr%s' % i
		numToChrom[i] = chrom
	for i, let in enumerate(letRange):
		chrom = 'chr%s' % let
		numToChrom[i + 1000] = chrom
	return numToChrom

def wigSort(fN):
	#sort by chromosome (123) and then XYZ...
	
	#init
	f = open(fN, 'r')
	header = f.readline() #header
	lines = f.readlines()
	f.close()
	chromToNum = returnChromToNumDict()
	numToChrom = returnNumToChromDict()
	
	#switch it to something sortable
	data = []
	for line in lines:
		data.append(line.split('\t'))
		data[-1][0] = chromToNum[data[-1][0]]
		
	#sort
	data = sorted(data, key=itemgetter(0), reverse = False) #sort by i
	
	#switch it back
	for i, line in enumerate(data):
		#print line
		data[i][0] = numToChrom[data[i][0]]
		data[i] = '\t'.join(data[i])
		
	
	f = open(fN, 'w')
	f.write(header)
	f.writelines(data)
	f.close()

def sortLines(sortLists, priorityList, convertList, orderList):
	'''Sortlists is a list of lists that will be sorted
	priority list is a list of ints with the column priorities (NOTE: Highest Priority is Last and Lowest is FIRST)
	convert list is a list of int that describes how each column should be converted before sorting
	orderList is a list of True/False values that describes if it should be top/bottom sorted
	Assumes that if you convert to a number it was originally a string.
	
	usage:
	pl = [0,1,3]
	cl = [1,1,1] # all ints
	ol = [True, True, True]
	
	sortLines(lines, pl, cl, ol)
	'''
	
	#convert string/ints: might have to keep track of what type they are...
	for i, val in enumerate(convertList): #different vals for int, float...
		if val == 0: continue
		
		column = priorityList[i]
		
		
		for l in sortLists:
			if val == 1: #int
				l[column] = int(l[column])
			elif val == 2:#float
				l[column] = float(l[column])
			elif val == 3: #string
				l[column] = str(l[column])
	
	#Now sort in correct priority.
	for i,item in enumerate(priorityList):
		#print 'sort step:'
		sortLists = sorted(sortLists, key=itemgetter(item), reverse = orderList[i])
		#print sortLists
	
	#Reconvert
	for i, val in enumerate(convertList):
		if val == 0: continue
		
		column = priorityList[i]
		
		for l in sortLists:
			if val == 1: #convert back to string
				l[column] = str(l[column])
			elif val == 2: #convert back to string
				l[column] = str(l[column])
			elif val == 3: #try int first and then float
				try:
					l[column] = int(l[column])
				except:
					l[column] = float(l[column])
				
	
	return sortLists

def sortLines(lists, column, colType = 'int', rev = False, toString = False):
		
	#establish original type
	originalType = type(lists[0][column])
	
	#establish conversions
	if colType == 'int':
		checkType = type(int())
		convertFxn = int
	
	#check if conversion needed and then convert	
	if originalType != checkType:
		for list in lists:
			list[column] = convertFxn(list[column])
	
	#sort it by column
	lists = sorted(lists, key=itemgetter(column), reverse = rev )

		
	#convert to string if need be
	if toString:
		for list in lists:
			list[column] = str(list[column])
			
	
	return lists
	
	
def sortFile(fN, column, type = 'int', rev = False):
	f = open(fN, 'r')
	lines = []
	for line in f:
		lines.append(line.strip().split('\t'))
	f.close()
	
	s = sortLines(lines, column, type, rev, toString = True)
	
	f = open(fN, 'w')
	for l in s:
		f.write('%s\n' % '\t'.join(l))
	f.close()

	

if __name__ == "__main__":
	import sys
	
	wigSort(sys.argv[1])
	
